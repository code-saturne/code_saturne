/*============================================================================
 * Checkpoint/restart handling for default application.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_physical_constants.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_restart_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_restart_default.c
        Checkpoint/restart handling for default application.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

const char *_coeff_name[] = {"bc_coeffs::a", "bc_coeffs::b",
                             "bc_coeffs::af", "bc_coeffs::bf",
                             "bc_coeffs::ad", "bc_coeffs::bd",
                             "bc_coeffs::ac", "bc_coeffs::bc"};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read and rebuild partial field metadata from legacy checkpoint.
 *
 * Note that when reading legacy files (Code_Saturne version 3.3 and below),
 * the old id will actually be the old scalar id (-1 for others).
 *
 * parameters:
 *   r <-- associated restart file pointer
 *----------------------------------------------------------------------------*/

static void
_read_legacy_field_info(cs_restart_t  *r)
{
  int retcode;

  int n_fields = cs_field_n_fields();

  /* Initialization */

  int kold = cs_field_key_id_try("old_scalar_num");

  /* Now read headers */

  cs_lnum_t n_old[4] = {0, 0, 0, 0}, n_cur[4] = {0, 0, 0, 0};

  const char *sec_id[] = {"nombre_variables",
                          "nombre_scalaires",
                          "nombre_scalaires_us",
                          "nombre_scalaires_pp"};

  for (int i = 0; i < 4; i++) {
    retcode = cs_restart_read_section(r,
                                      sec_id[i],
                                      CS_MESH_LOCATION_NONE,
                                      1,
                                      CS_TYPE_cs_int_t,
                                      n_old + i);
    if (retcode != CS_RESTART_SUCCESS)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Error reading variable information in restart file \"%s\"."),
         cs_restart_get_name(r));
  }

  const int kv = cs_field_key_id_try("variable_id");
  const int ks = cs_field_key_id_try("scalar_id");

  /* Count variables and user and non-user scalars */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *f = cs_field_by_id(f_id);
    int v_num = -1;
    if (kv > -1)
      v_num = cs_field_get_key_int(f, kv);
    if (v_num > 0) {
      int s_num = -1;
      n_cur[0] += 1;
      if (ks > -1)
        s_num = cs_field_get_key_int(f, ks);
      if (s_num > -1) {
        n_cur[1] += 1;
        if (f->type & CS_FIELD_USER)
          n_cur[2] += 1;
        else
          n_cur[3] += 1;
      }
    }
  }

  /* Possible shift in old ids if temperature has been moved from
     user to model scalar */

  int us_shift = 0, pp_shift = 0;

  /* Special case if temperature has been moved from
     user to model scalar */

  if (   n_cur[1] == n_old[1] && n_cur[2] == n_old[2] -1
      && n_cur[3] == 1 && n_old[3] == 0) {

    if (CS_F_(t) != NULL || CS_F_(h) != NULL) {
      us_shift = -1;
      pp_shift = n_cur[2];
    }

  }

  /* Warn in case of change */

  if (   n_cur[0] != n_old[0] || n_cur[1] != n_old[1]
      || n_cur[2] != n_old[2] || n_cur[3] != n_old[3]) {

    /* Special case if temperature has been moved from
       user to model scalar */

    if (n_cur[0] == n_old[0] && n_cur[1] == n_old[1] && us_shift == -1)
      bft_printf
        (_("\nRemark: the thermal scalar was treated as a user scalar\n"
           "          in the restart file, and is moved to a model scalar\n"
           "          in the current computation.\n"));

    else {
      bft_printf
        (_("\n"
           "  Warning: the number of variables or scalars has been changed\n"
           "           relative to the restart file.\n\n"
           "  currently  %d variables, of which %d scalars\n"
           "  previously %d variables, of which %d scalars\n\n"
           "  The computation continues, with a partial restart.\n"),
         (int)(n_cur[0]), (int)(n_cur[1]), (int)(n_old[0]), (int)(n_old[1]));

    }

  }

  /* Now check (and update if necessary) old scalar id */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int old_scal_num = -1;
    int s_num = -1;
    if (ks > -1)
      s_num = cs_field_get_key_int(f, ks);
    if (s_num > -1) {
      old_scal_num = -1;
      if (kold > -1)
        old_scal_num = cs_field_get_key_int(f, kold);
      if (old_scal_num < 0) {
        if (f->type & CS_FIELD_USER)
          old_scal_num = s_num + us_shift;
        else
          old_scal_num = s_num + pp_shift;
        if (old_scal_num > n_old[1])
          old_scal_num = -1;
      }
      else {
        if (old_scal_num > n_old[1])
          bft_error
            (__FILE__, __LINE__, 0,
             _("Field \"%s\" has user-defined key \"old_scalar_num\" value %d,\n"
               "but the number of available scalars in restart is %d."),
             f->name, old_scal_num, (int)(n_old[1]));
      }
      if (kold < 0)
        kold = cs_field_define_key_int("old_scalar_num",
                                       -1,
                                       CS_FIELD_VARIABLE);
      cs_field_set_key_int(f, kold, old_scal_num);
    }
  }

}

/*----------------------------------------------------------------------------
 * Synchronize cell-based field values.
 *
 * parameters:
 *   f    <-> field whose values should be synchronized
 *   t_id <-- time id (0 for current, 1 for previous, ...)
 *----------------------------------------------------------------------------*/

static void
_sync_field_vals(cs_field_t  *f,
                 int          t_id)
{
  const cs_mesh_t *m = cs_glob_mesh;

  if (m->halo != NULL) {

    cs_halo_type_t  halo_type = CS_HALO_EXTENDED;
    cs_real_t      *v = f->vals[t_id];

    cs_halo_sync_var_strided(m->halo, halo_type, v, f->dim);

    if (m->n_init_perio > 0) {
      if (f->dim == 3)
        cs_halo_perio_sync_var_vect(m->halo, halo_type, v, 3);
      else if (f->dim == 6)
        cs_halo_perio_sync_var_sym_tens(m->halo, halo_type, v);
      else if (f->dim == 9)
        cs_halo_perio_sync_var_tens(m->halo, halo_type, v);
      else if (f->dim == 1 && f == CS_F_(r13)) {
        cs_halo_perio_sync_var_tens_ni(m->halo,
                                       halo_type,
                                       CS_F_(r11)->vals[t_id],
                                       CS_F_(r12)->vals[t_id],
                                       CS_F_(r13)->vals[t_id],
                                       CS_F_(r12)->vals[t_id],
                                       CS_F_(r22)->vals[t_id],
                                       CS_F_(r23)->vals[t_id],
                                       CS_F_(r13)->vals[t_id],
                                       CS_F_(r23)->vals[t_id],
                                       CS_F_(r33)->vals[t_id]);
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Read field values from checkpoint.
 *
 * Values are found using the default rules based on the field's name
 * postfixed by ::vals::%t_id, or its name itself.
 *
 * parameters:
 *   r            <-- associated restart file pointer
 *   restart_name <-- base name with which read is attempted
 *   t_id         <-- time id (0 for current, 1 for previous, ...)
 *   f            <-> field whose values should be read
 *
 * returns:
 *   CS_RESTART_SUCCESS in case of success, CS_RESTART_ERR_... otherwise
 *----------------------------------------------------------------------------*/

static int
_read_field_vals(cs_restart_t  *r,
                 const char    *r_name,
                 int            t_id,
                 cs_field_t    *f)
{
  int retcode = CS_RESTART_SUCCESS;

  char _sec_name[128];
  char *sec_name = _sec_name;

  if (strlen(r_name) > 96)
    BFT_MALLOC(sec_name, strlen(r_name) + 64, char); /* wide margin */

  /* Check for data; data will be read later, so that compatibility
     checks may be done first; we really try reading the data only
     at the end, so if it is not found, a warning will be logged only
     once (and not once per test), and preferentially use the
     base (non-compatibility) name. */

  snprintf(sec_name, 127, "%s::vals::%d", r_name, t_id);
  sec_name[127] = '\0';

  retcode = cs_restart_check_section(r,
                                     sec_name,
                                     f->location_id,
                                     f->dim,
                                     CS_TYPE_cs_real_t);

  /* Otherwise, try reading with basic (restart) name only if requested */

  if (   (retcode == CS_RESTART_ERR_EXISTS || retcode == CS_RESTART_ERR_N_VALS)
      && r_name != f->name) {
    snprintf(sec_name, 127, "%s", r_name);
    sec_name[127] = '\0';
    retcode = cs_restart_check_section(r,
                                       sec_name,
                                       f->location_id,
                                       f->dim,
                                       CS_TYPE_cs_real_t);
  }

  /* Read if available */

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(r,
                                      sec_name,
                                      f->location_id,
                                      f->dim,
                                      CS_TYPE_cs_real_t,
                                      f->vals[t_id]);

  /* Try to read anyways (with base name) to log warning */

  else {
    snprintf(sec_name, 127, "%s::vals::%d", r_name, t_id);
    sec_name[127] = '\0';
    retcode = cs_restart_read_section(r,
                                      sec_name,
                                      f->location_id,
                                      f->dim,
                                      CS_TYPE_cs_real_t,
                                      f->vals[t_id]);
  }

  if (sec_name != _sec_name)
    BFT_FREE(sec_name);

  if (   retcode == CS_RESTART_SUCCESS
      && f->location_id == CS_MESH_LOCATION_CELLS)
    _sync_field_vals(f, t_id);

  return retcode;
}

/*----------------------------------------------------------------------------
 * Read field values from legacy checkpoint.
 *
 * Values are found using older names for compatibility with older files.
 * For cell-based fields, the old name base is appended automatically with
 * "_ce_phase01", except for scalars, where the name uses a different scheme,
 * based on "scalaire_ce_%04" % s_num.
 *
 * parameters:
 *   r            <-- associated restart file pointer
 *   r_name <-- base name with which read is attempted
 *   t_id         <-- time id (0 for current, 1 for previous, ...)
 *   f            <-> file whose values should be read
 *
 * returns:
 *   CS_RESTART_SUCCESS in case of success, CS_RESTART_ERR_... otherwise
 *----------------------------------------------------------------------------*/

static int
_read_field_vals_legacy(cs_restart_t  *r,
                        const char    *r_name,
                        int            t_id,
                        cs_field_t    *f)
{
  char sec_name[128] = "";
  char old_name_x[128] = "", old_name_y[128] = "", old_name_z[128] = "";
  char old_name_xx[128] = "", old_name_yy[128] = "", old_name_zz[128] = "";
  char old_name_xy[128] = "", old_name_yz[128] = "", old_name_xz[128] = "";

  int retcode = CS_RESTART_SUCCESS;

  /* Check for renaming */

  char old_name[128] = "";
  int ks = cs_field_key_id_try("scalar_id");
  int scalar_id = cs_field_get_key_int(f, ks);

  /* Special case for scalars */

  if (scalar_id > -1) {
    if (r_name != f->name) {
      const char *name = r_name;
      while (*name != '\0' && !isdigit(*name))
        name++;
      scalar_id = atoi(name) - 1;
    }
    if (scalar_id > -1)
      snprintf(old_name, 127, "%04d", scalar_id+1);
    else
      snprintf(old_name, 127, "%s", r_name);
  }

  /* Other fields may need specific renaming
     (old_name for partial section name, sec_name for direct section name) */

  else if (r_name == f->name) {
    snprintf(old_name, 127, "%s", f->name);
    if (f == CS_F_(vel)) {
      if (t_id == 0)
        strncpy(old_name, "vitesse", 127);
      else if (t_id == 1)
        strncpy(sec_name, "velocity_prev", 127);
    }
    else if (f == CS_F_(p))
      strncpy(old_name, "pression", 127);
    else if (f == CS_F_(r11))
      strncpy(old_name, "R11", 127);
    else if (f == CS_F_(r22))
      strncpy(old_name, "R22", 127);
    else if (f == CS_F_(r33))
      strncpy(old_name, "R33", 127);
    else if (f == CS_F_(r12))
      strncpy(old_name, "R12", 127);
    else if (f == CS_F_(r13))
      strncpy(old_name, "R13", 127);
    else if (f == CS_F_(r23))
      strncpy(old_name, "R23", 127);
    else if (f == CS_F_(rij))
      strncpy(old_name, "Rij", 127);
    else if (f == CS_F_(eps))
      strncpy(old_name, "eps", 127);
    else if (f == CS_F_(f_bar))
      strncpy(old_name, "fb", 127);
    else if (f == CS_F_(alp_bl)) {
      /* Special case: "al" also possible here, depending on turbulence model;
         check for either, with one test for main restart, the other for
         the auxilairy restart */
      int sec_code;
      strncpy(old_name, "alp", 127);
      sec_code = cs_restart_check_section(r, "al_ce_phase01",
                                          1, 1, CS_TYPE_cs_real_t);
      if (sec_code == CS_RESTART_SUCCESS)
        strncpy(old_name, "al", 127);
      else
        sec_code = cs_restart_check_section(r, "fm_al_phase01",
                                            0, 1, CS_TYPE_cs_int_t);
      if (sec_code == CS_RESTART_SUCCESS)
        strncpy(old_name, "al", 127);
    }
    else if (f == CS_F_(nusa))
      strncpy(old_name, "nusa", 127);
    else if (f == CS_F_(mesh_u))
      strncpy(old_name, "vit_maillage", 127);
    else if (f == CS_F_(rho)) {
      if (t_id == 0)
        strncpy(old_name, "rho", 127);
      else if (t_id == 1)
        strncpy(old_name, "rho_old", 127);
    }
    else if (f == CS_F_(rho_b))
      strncpy(sec_name, "rho_fb_phase01", 127);

    else if (f == CS_F_(cp))
      strncpy(old_name, "cp", 127);

    else if (f == CS_F_(mu))
      strncpy(old_name, "viscl", 127);
    else if (f == CS_F_(mu_t))
      strncpy(old_name, "visct", 127);

    else if (f == CS_F_(t_b))
      strncpy(old_name, "tparoi_fb", 127);
    else if (f == CS_F_(qinci))
      strncpy(old_name, "qincid_fb", 127);
    else if (f == CS_F_(hconv))
      strncpy(old_name, "hfconv_fb", 127);
    else if (f == CS_F_(fconv))
      strncpy(old_name, "flconv_fb", 127);

    else if (strcmp(f->name, "dt") == 0)
      strncpy(sec_name, "dt_variable_espace_ce", 127);

    else if (strcmp(f->name, "volume_forces") == 0)
      strncpy(sec_name, "force_ext_ce_phase01", 127);

    else if (strcmp(f->name, "hydrostatic_pressure_prd") == 0)
      strncpy(sec_name, "Prhyd_pre_phase01", 127);

    else if (strcmp(f->name, "void_fraction") == 0)
      strncpy(sec_name, "taux_vide_ce", 127);

    else if (strcmp(f->name, "rad_st") == 0)
      strncpy(sec_name, "rayexp_ce", 127);
    else if (strcmp(f->name, "rad_st_implicit") == 0)
      strncpy(sec_name, "rayimp_ce", 127);
    else if (f == CS_F_(rad_lumin))
      strncpy(sec_name, "luminance", 127);

    else if (strcmp(f->name, "joule_power") == 0)
      strncpy(sec_name, "tsource_sc_ce_joule", 127);
    else if (strcmp(f->name, "laplace_force") == 0)
      strncpy(old_name, "laplace_force", 127);
  }

  if (sec_name[0] == '\0') {
    if (scalar_id > -1)
      snprintf(sec_name, 127, "scalaire_ce_%04d", scalar_id);
    else if (f->location_id == CS_MESH_LOCATION_CELLS)
      snprintf(sec_name, 127, "%s_ce_phase01", old_name);
    else
      snprintf(sec_name, 127, "%s", old_name);
  }

  retcode = cs_restart_check_section(r,
                                     sec_name,
                                     f->location_id,
                                     f->dim,
                                     CS_TYPE_cs_real_t);

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(r,
                                      sec_name,
                                      f->location_id,
                                      f->dim,
                                      CS_TYPE_cs_real_t,
                                      f->vals[t_id]);

  /* Last chance for 3D fields */

  else if (f->dim == 3 && retcode == CS_RESTART_ERR_EXISTS) {

    if (strcmp(old_name, "vit_maillage") == 0) {
      snprintf(old_name_x, 127, "%s_u_ce", old_name);
      snprintf(old_name_y, 127, "%s_v_ce", old_name);
      snprintf(old_name_z, 127, "%s_w_ce", old_name);
    }
    else if (strcmp(old_name, "laplace_force") == 0) {
      snprintf(old_name_x, 127, "%s_1", old_name);
      snprintf(old_name_y, 127, "%s_2", old_name);
      snprintf(old_name_z, 127, "%s_2", old_name);
    }
    else {
      snprintf(old_name_x, 127, "%s_u_ce_phase01", old_name);
      snprintf(old_name_y, 127, "%s_v_ce_phase01", old_name);
      snprintf(old_name_z, 127, "%s_w_ce_phase01", old_name);
    }

    retcode = cs_restart_check_section(r,
                                       old_name_x,
                                       f->location_id,
                                       1,
                                       CS_TYPE_cs_real_t);

    if (retcode == CS_RESTART_SUCCESS)
      retcode = cs_restart_read_real_3_t_compat(r,
                                                sec_name,
                                                old_name_x,
                                                old_name_y,
                                                old_name_z,
                                                f->location_id,
                                                (cs_real_3_t *)(f->vals[t_id]));
  }
  else if (f->dim == 6 && retcode == CS_RESTART_ERR_EXISTS) {

    if (strcmp(old_name, "Rij") == 0) {
      snprintf(old_name_xx, 127, "r11::vals::0");
      snprintf(old_name_yy, 127, "r22::vals::0");
      snprintf(old_name_zz, 127, "r33::vals::0");
      snprintf(old_name_xy, 127, "r12::vals::0");
      snprintf(old_name_yz, 127, "r23::vals::0");
      snprintf(old_name_xz, 127, "r13::vals::0");

      retcode = cs_restart_check_section(r,
                                         old_name_xx,
                                         f->location_id,
                                         1,
                                         CS_TYPE_cs_real_t);

      if (retcode == CS_RESTART_SUCCESS)
        retcode = cs_restart_read_real_6_t_compat(r,
                                                  "rij::vals::0",
                                                  old_name_xx,
                                                  old_name_yy,
                                                  old_name_zz,
                                                  old_name_xy,
                                                  old_name_yz,
                                                  old_name_xz,
                                                  f->location_id,
                                                  (cs_real_6_t *)(f->vals[t_id]));
    }
  }

  if (   retcode == CS_RESTART_SUCCESS
      && f->location_id == CS_MESH_LOCATION_CELLS)
    _sync_field_vals(f, t_id);

  return retcode;
}

/*----------------------------------------------------------------------------
 * Determine mass flux number of a given variable for legacy restart file
 *
 * parameters:
 *   r          <-> associated restart file pointer
 *   f          <-- associated field pointer
 *   scalar_num <-- associated scalar number, or -1
 *   t_id       <-- associated time id (0: current, 1: previous)
 *
 * returns:
 *   number of matching mass flux in restart file (-1 if read failed)
 *----------------------------------------------------------------------------*/

static int
_legacy_mass_flux_num(cs_restart_t      *r,
                      const cs_field_t  *f,
                      int                scalar_num,
                      int                t_id)
{
  int retval = 1;

  /* As of Code_Saturne 3.3, only scalars may have a different mass flux
     from the "main" mass flux (in the case of scalars with drift), so for
     all others, reading the associated mass flux name is of no real use. */

  char sec_name[128] = "";

  const char *prefix[2] = {"fm_", "fm_a_"};
  if (scalar_num > 0)
    snprintf(sec_name, 127, "%sscalaire%04d", prefix[t_id], scalar_num);
  else if (strcmp(f->name, "void_fraction") == 0)
    snprintf(sec_name, 127, "%staux_vide", prefix[t_id]);

  /* Read from restart */

  if (sec_name[0] != '\0') {
    cs_lnum_t buf[1];
    sec_name[127] = '\0';
    int retcode = cs_restart_read_section(r,
                                          sec_name,
                                          CS_MESH_LOCATION_NONE,
                                          1,
                                          CS_TYPE_cs_int_t,
                                          buf);
    if (retcode == CS_RESTART_SUCCESS)
      retval = buf[0];
    else
      retval = -1;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Read fields depending on others from checkpoint.
 *
 * This function handles legacy files (Code_Saturne version 3.3 and below).
 *
 * parameters:
 *   r         <-> associated restart file pointer
 *   key       <-> key for field association
 *   read_flag <-> flag to track fields read, or NULL;
 *                 set to 1 for fields read, -1 for fields
 *                 failed to read (size: n_fields)
 *
 * returns:
 *   number of fields read
 *----------------------------------------------------------------------------*/

static int
_read_linked_fields_legacy(cs_restart_t  *r,
                           const char    *key,
                           int            read_flag[])
{
  int retcode;

  int retcount = 0;

  /* Initialization */

  int category = 0;

  const int n_fields = cs_field_n_fields();

  const int key_id = cs_field_key_id(key);
  const int key_flag = cs_field_key_flag(key_id);

  const int kold = cs_field_key_id_try("old_scalar_num");
  const int ks = cs_field_key_id_try("scalar_id");

  /* Determine field type (out of possibilities in legacy files) */

  if (strcmp(key, "inner_mass_flux_id") == 0)
    category = 1;
  else if (strcmp(key, "boundary_mass_flux_id") == 0)
    category = 2;
  else if (strcmp(key, "diffusivity_id") == 0)
    category = 3;

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t *f = cs_field_by_id(f_id);

    if (key_flag == -1 || !(f->type & key_flag))
      continue;

    const int lnk_f_id = cs_field_get_key_int(f, key_id);
    int s_num = -1;

    if (lnk_f_id > -1) {

      cs_field_t *f_lnk = cs_field_by_id(lnk_f_id);

      if (read_flag[lnk_f_id] != 0)
        continue;

      read_flag[lnk_f_id] = -1;

      /* check for (possibly renumbered) scalar */
      if (f->type & CS_FIELD_VARIABLE) {
        if (kold > -1)
          s_num = cs_field_get_key_int(f, kold);
        if (s_num < 0 && ks > -1)
          s_num = cs_field_get_key_int(f, ks);
      }

      for (int t_id = 0; t_id < 2; t_id++) {

        if (t_id >= f_lnk->n_time_vals)
          break;

        /* Build field name to read based on category */

        char sec_name[128];
        if (category == 1) {
          int mf_num = _legacy_mass_flux_num(r, f, s_num, t_id);
          if (t_id == 0)
            snprintf(sec_name, 127, "flux_masse_fi_%04d", mf_num);
          else
            snprintf(sec_name, 127, "flux_masse_a_fi_%04d", mf_num);
        }
        else if (category == 2) {
          int mf_num = _legacy_mass_flux_num(r, f, s_num, t_id);
          if (t_id == 0)
            snprintf(sec_name, 127, "flux_masse_fb_%04d", mf_num);
          else
            snprintf(sec_name, 127, "flux_masse_a_fb_%04d", mf_num);
        }
        else if (category == 3)
          snprintf(sec_name, 127, "visls_ce_scalaire%04d", s_num);

        /* Now we know which field name to read */

        retcode = cs_restart_check_section(r,
                                           sec_name,
                                           f->location_id,
                                           f->dim,
                                           CS_TYPE_cs_real_t);

        if (retcode == CS_RESTART_SUCCESS)
          retcode = cs_restart_read_section(r,
                                            sec_name,
                                            f->location_id,
                                            f->dim,
                                            CS_TYPE_cs_real_t,
                                            f->vals[t_id]);

        if (retcode == CS_RESTART_SUCCESS) {
          if (t_id == 0)
            read_flag[lnk_f_id] = 1;
          else
            read_flag[lnk_f_id] += 2;
          retcount += 1;
        }

      } /* t_id */

    }

  }

  return retcount;
}

/*----------------------------------------------------------------------------
 * Compare old and new values for a given field key
 *
 * parameters:
 *   r             <-- associated restart file pointer
 *   key           <-- associated model key
 *   old_field_map <-- name to id map of fields in restart file
 *
 * returns:
 *   number of values changed, or -1 if information was not found
 *----------------------------------------------------------------------------*/

static int
_check_field_model(cs_restart_t               *r,
                   const char                 *key,
                   const cs_map_name_to_id_t  *old_field_map)
{
  int retcode;

  const int n_fields = cs_field_n_fields();
  const int n_o_fields = cs_map_name_to_id_size(old_field_map);

  const int key_id = cs_field_key_id(key);
  const int key_flag = cs_field_key_flag(key_id);
  const int kr = cs_field_key_id_try("restart_name");

  int  n_diff = 0;

  cs_lnum_t *old_key_val;
  BFT_MALLOC(old_key_val, n_o_fields, cs_lnum_t);

  char *sec_name;
  BFT_MALLOC(sec_name, strlen("fields:") + strlen(key) + 1, char);
  strcpy(sec_name, "fields:");
  strcat(sec_name, key);

  /* Read metadata */

  retcode = cs_restart_check_section(r,
                                     sec_name,
                                     CS_MESH_LOCATION_NONE,
                                     n_o_fields,
                                     CS_TYPE_cs_int_t);

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(r,
                                      sec_name,
                                      CS_MESH_LOCATION_NONE,
                                      n_o_fields,
                                      CS_TYPE_cs_int_t,
                                      old_key_val);

  /* If data is available, compare models */

  if (retcode == CS_RESTART_SUCCESS) {

    for (int f_id = 0; f_id < n_fields; f_id++) {

      const cs_field_t *f = cs_field_by_id(f_id);

      if (key_flag == -1 || !(f->type & key_flag))
        continue;

      const char *f_name = NULL;
      if (kr > -1)
        f_name = cs_field_get_key_str(f, kr);
      if (f_name == NULL)
        f_name = f->name;

      int old_f_id = cs_map_name_to_id_try(old_field_map, f_name);
      if (old_f_id > -1) {
        if (cs_field_get_key_int(f, key_id) != old_key_val[old_f_id])
          n_diff += 1;
      }

    }

  }
  else if (retcode == CS_RESTART_ERR_EXISTS)
    n_diff = -1;

  else
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error %d reading \"%s\" in restart file \"%s\"."),
       retcode, sec_name, cs_restart_get_name(r));

  BFT_FREE(sec_name);
  BFT_FREE(old_key_val);

  return n_diff;
}

/*----------------------------------------------------------------------------
 * Check old turbulent flux model
 *
 * parameters:
 *   r             <-- associated restart file pointer
 *   old_field_map <-- name to id map of fields in restart file
 *----------------------------------------------------------------------------*/

static void
_check_turb_flux_model(cs_restart_t               *r,
                       const cs_map_name_to_id_t  *old_field_map)
{
  const int n_fields = cs_field_n_fields();

  int  n_diff = _check_field_model(r,
                                   "turbulent_flux_model",
                                   old_field_map);

  /* Read in legacy mode if required */

  if (n_diff < 0) {

    int kold = cs_field_key_id_try("old_scalar_num");
    const int key_id = cs_field_key_id("turbulent_flux_model");

    if (kold > -1) {

      n_diff = 0;

      for (int f_id = 0; f_id < n_fields; f_id++) {

        const cs_field_t *f = cs_field_by_id(f_id);

        if (!(f->type & CS_FIELD_VARIABLE))
          continue;

        int s_num = cs_field_get_key_int(f, kold);

        if (s_num > 0) {
          char sec_name[128];
          cs_lnum_t old_s_model[1];
          snprintf(sec_name, 127, "turbulent_flux_model%04d", s_num);
          sec_name[127] = '\0';
          int retcode = cs_restart_read_section(r,
                                                sec_name,
                                                CS_MESH_LOCATION_NONE,
                                                1,
                                                CS_TYPE_cs_int_t,
                                                old_s_model);
          if (retcode == CS_RESTART_SUCCESS) {
            if (cs_field_get_key_int(f, key_id) != old_s_model[0])
              n_diff += 1;
          }
        }

      }

    }

  }

  if (n_diff > 0)
    bft_printf
      (_("\n"
         "  Warning: the turbulent flux model has been changed\n"
         "           for %d fields relative to the restart file\n\n"
         "  The computation continues, with a partial restart.\n"),
       n_diff);
}

/*----------------------------------------------------------------------------
 * Read model option from file with compatibility for older files.
 *
 * parameters:
 *   r          <-- associated restart file pointer
 *   m_name     <-- name of model in current vesions of the restart file
 *   m_name_old <-- name of model in different vesions of the restart file
 *   count      <-- number of values with current name
 *   count_old  <-- number of values with older name
 *   options    <-> model flags
 *
 * returns:
 *   number of values read
 *----------------------------------------------------------------------------*/

static int
_read_model_option_compat(cs_restart_t  *r,
                          const char    *m_name,
                          const char    *m_name_old,
                          int            count,
                          int            count_old,
                          cs_lnum_t      options[])
{
  int retcount = 0;

  int retcode = cs_restart_check_section(r,
                                         m_name,
                                         CS_MESH_LOCATION_NONE,
                                         1,
                                         CS_TYPE_cs_int_t);

  if (retcode == CS_RESTART_ERR_EXISTS) {
    retcode = cs_restart_check_section(r,
                                       m_name_old,
                                       CS_MESH_LOCATION_NONE,
                                       count_old,
                                       CS_TYPE_cs_int_t);
    if (retcode == CS_RESTART_SUCCESS) {
      cs_restart_read_section(r,
                              m_name_old,
                              CS_MESH_LOCATION_NONE,
                              count_old,
                              CS_TYPE_cs_int_t,
                              options);
      if (retcode == CS_RESTART_SUCCESS)
        retcount = count_old;
    }
  }
  else {
    cs_restart_read_section(r,
                            m_name,
                            CS_MESH_LOCATION_NONE,
                            count,
                            CS_TYPE_cs_int_t,
                            options);
    if (retcode == CS_RESTART_SUCCESS)
      retcount = count;
  }

  return retcount;
}

/*----------------------------------------------------------------------------
 * Read turbulence variable section to a 1d array of cells.
 *
 * This wrapper function is intended mainly to allow calls with simplifed
 * arguments.
 *
 * parameters:
 *   r         <-- associated restart file pointer
 *   base_name <-- base name of section to read
 *   old_name  <-- old name of section to read
 *   t_id      <-- associated time id
 *   vals      --> pointer to values
 *
 * returns:
 *   1 in case of error, 0 in case of success
 *----------------------------------------------------------------------------*/

static int
_read_turb_array_1d_compat(cs_restart_t   *r,
                           const char     *base_name,
                           const char     *old_name,
                           int             t_id,
                           cs_real_t      *vals)
{
  char sec_name[128];
  snprintf(sec_name, 127, "%s::vals::%01d", base_name, t_id);
  sec_name[127] = '\0';
  int retcode = CS_RESTART_SUCCESS;
  int retval = 0;

  if (t_id == 0) {
    char old_sec_name[128];
    snprintf(old_sec_name, 127, "%s_ce_phase01", old_name);
    old_sec_name[127] = '\0';
    retcode = cs_restart_read_section_compat(r,
                                             sec_name,
                                             old_sec_name,
                                             CS_MESH_LOCATION_CELLS,
                                             1,
                                             CS_TYPE_cs_real_t,
                                             vals);
  }
  else
    retcode = cs_restart_read_section(r,
                                      sec_name,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      CS_TYPE_cs_real_t,
                                      vals);
  if (retcode != CS_RESTART_SUCCESS)
    retval = 0;

  return retval;
}

/*----------------------------------------------------------------------------
 * Compare current and previous turbulence model.
 *
 * parameters:
 *   r         <-- associated restart file pointer
 *
 * returns:
 *   previous turbulence model
 *----------------------------------------------------------------------------*/

static int
_read_turbulence_model(cs_restart_t  *r)
{
  cs_lnum_t iturb_old = 0;

  /* Determine previous turbulence model */

  int n_opt_vals = _read_model_option_compat(r,
                                             "turbulence_model",
                                             "modele_turbulence_phase01",
                                             1,
                                             1,
                                             &iturb_old);

  if (n_opt_vals != 1)
    iturb_old = -999;

  return iturb_old;
}

/*----------------------------------------------------------------------------
 * Read turbulence variables converted from another model.
 *
 * This function should be called once we have read all "matching"
 * variables, as it assumes variables common to several turbulence models
 * (for example k between k-epsilon and k-omega, or epsilon between
 * k-epsilon and Rij-epsilon) have already been read.
 *
 * parameters:
 *   r         <-- associated restart file pointer
 *   iturb_cur <-- current turbulence model number
 *   iturb_old <-- previous turbulence model number
 *   t_id      <-- associated time id
 *   read_flag <-> flag to track fields read, or NULL;
 *                 set to 1 for fields read, -1 for fields
 *                 failed to read (size: n_fields)
 *----------------------------------------------------------------------------*/

static void
_read_and_convert_turb_variables(cs_restart_t  *r,
                                 int            iturb_cur,
                                 int            iturb_old,
                                 int            t_id,
                                 int            read_flag[])
{
  int   err_sum = 0;

  const int itytur_cur = iturb_cur/10;
  const int itytur_old = iturb_old/10;

  /* If the turbulence modela has not changed, nothing to do */

  if (iturb_cur == iturb_old)
    return;

  /* Warn user that turbulence model has changed */

  bft_printf
    (_("\n"
       "  Warning: the turbulence model has been changed\n"
       "           relative to the restart file.\n\n"
       "  current model:  %d\n"
       "  previous model: %d\n\n"
       "  The computation continues, with a partial an/or adapted restart.\n"),
     iturb_cur, (int)iturb_old);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_real_t *v_tmp;
  BFT_MALLOC(v_tmp, n_cells, cs_real_t);

  /* Now convert variables if needed. When we do not know how to
     deduce turbulence variables from available info, we ignore
     them, and keep the default initializations. */

  int t_mask = (t_id == 0) ? 1 : 2 << (t_id-1);

  if (itytur_cur == 2) { /* k-eps or v2f (phi-fbar or BL-v2/k) to k-eps */

    /* if (itytur_old == 5),
       k and epsilon are common (so already read), and phi or f-bar
       is ignored, so nothing left to do here */

    if (itytur_old == 3) { /* Rij to k-epsilon (epsilon already read) */

      cs_real_t *v_k = CS_F_(k)->vals[t_id];

      err_sum += _read_turb_array_1d_compat(r, "r11", "R11", t_id, v_k);

      err_sum += _read_turb_array_1d_compat(r, "r22", "R22", t_id, v_tmp);

      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_k[i] += v_tmp[i];

      err_sum += _read_turb_array_1d_compat(r, "r33", "R33", t_id, v_tmp);

      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_k[i] = 0.5 * (v_k[i] + v_tmp[i]);

      if (err_sum == 0)
        read_flag[CS_F_(k)->id] += t_mask;

    } else if (iturb_old == 60) { /* k-omega to k-epsilon (k already read) */

      cs_real_t *v_eps = CS_F_(eps)->vals[t_id];
      const cs_real_t *v_k = CS_F_(eps)->vals[t_id];

      err_sum += _read_turb_array_1d_compat(r, "omega", "omega", t_id, v_eps);

      /* Now transform omega to epsilon */

      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_eps[i] = cs_turb_cmu*v_eps[i]*v_k[i];

      if (err_sum == 0)
        read_flag[CS_F_(eps)->id] += t_mask;
    }

  } else if (itytur_cur == 3 && CS_F_(rij) != NULL) { /* New computation
                                                         is in Rij-epsilon */

    if (itytur_old == 2 || iturb_old == 50) { /* k-epsilon or v2f
                                                 (phi-fbar or BL-v2/k) to-> Rij
                                                 (epsilon already read) */

      cs_real_6_t *v_rij = (cs_real_6_t *)(CS_F_(rij)->vals[t_id]);
      cs_real_t *v_k;
      BFT_MALLOC(v_k, n_cells, cs_real_t);

      err_sum += _read_turb_array_1d_compat(r, "k", "k", t_id, v_k);

      double d2s3 = 2./3.;

      for (cs_lnum_t i = 0; i < n_cells; i++) {
        double d2s3xk = v_k[i]*d2s3;
        v_rij[i][0] = d2s3xk;
        v_rij[i][1] = d2s3xk;
        v_rij[i][2] = d2s3xk;
        v_rij[i][3] = 0.;
        v_rij[i][4] = 0.;
        v_rij[i][5] = 0.;
      }

      if (err_sum == 0)
        read_flag[CS_F_(rij)->id] += t_mask;

      BFT_FREE(v_k);

    }

    /* if (itytur_old == 3) Rij to Rij variant; nothing to do
     *                      (when switching to EBRSM to another model,
     *                      keep alpha to default initialization) */

    else if (iturb_old == 60) { /* k-omega to Rij */

      cs_real_6_t *v_rij = (cs_real_6_t *)(CS_F_(rij)->vals[t_id]);
      cs_real_t *v_eps = CS_F_(eps)->vals[t_id];

      cs_real_t *v_k;
      BFT_MALLOC(v_k, n_cells, cs_real_t);

      err_sum += _read_turb_array_1d_compat(r, "k", "k", t_id, v_k);
      err_sum += _read_turb_array_1d_compat(r, "omega", "omega", t_id, v_eps);

      /* Now transform omega to epsilon */

      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_eps[i] = cs_turb_cmu*v_eps[i]*v_k[i];

      double d2s3 = 2./3.;

      for (cs_lnum_t i = 0; i < n_cells; i++) {
        double d2s3xk = v_k[i]*d2s3;
        v_rij[i][0] = d2s3xk;
        v_rij[i][1] = d2s3xk;
        v_rij[i][2] = d2s3xk;
        v_rij[i][3] = 0.;
        v_rij[i][4] = 0.;
        v_rij[i][5] = 0.;
      }

      if (err_sum == 0)
        read_flag[CS_F_(rij)->id] += t_mask;

      BFT_FREE(v_k);
    }

  } else if (itytur_cur == 3) { /* New computation is in Rij-epsilon */
                                /* using the older, non-interleaved variant */

    if (itytur_old == 2 || iturb_old == 50) { /* k-epsilon or v2f
                                                 (phi-fbar or BL-v2/k) to-> Rij
                                                 (epsilon already read) */

      cs_real_t *v_r11 = CS_F_(r11)->vals[t_id];
      cs_real_t *v_r22 = CS_F_(r22)->vals[t_id];
      cs_real_t *v_r33 = CS_F_(r33)->vals[t_id];
      cs_real_t *v_r12 = CS_F_(r12)->vals[t_id];
      cs_real_t *v_r23 = CS_F_(r23)->vals[t_id];
      cs_real_t *v_r13 = CS_F_(r13)->vals[t_id];

      err_sum += _read_turb_array_1d_compat(r, "k", "k", t_id, v_r11);

      double d2s3 = 2./3.;

      for (cs_lnum_t i = 0; i < n_cells; i++) {
        double d2s3xk = v_r11[i]*d2s3;
        v_r11[i] = d2s3xk;
        v_r22[i] = d2s3xk;
        v_r33[i] = d2s3xk;
        v_r12[i] = 0.;
        v_r23[i] = 0.;
        v_r13[i] = 0.;
      }

      if (err_sum == 0) {
        read_flag[CS_F_(r11)->id] += t_mask;
        read_flag[CS_F_(r22)->id] += t_mask;
        read_flag[CS_F_(r33)->id] += t_mask;
        read_flag[CS_F_(r12)->id] += t_mask;
        read_flag[CS_F_(r23)->id] += t_mask;
        read_flag[CS_F_(r13)->id] += t_mask;
      }

    }

    /* if (itytur_old == 3) Rij to Rij variant; nothing to do
     *                      (when switching to EBRSM to another model,
     *                      keep alpha to default initialization) */

    else if (iturb_old == 60) { /* k-omega to Rij */

      cs_real_t *v_r11 = CS_F_(r11)->vals[t_id];
      cs_real_t *v_r22 = CS_F_(r22)->vals[t_id];
      cs_real_t *v_r33 = CS_F_(r33)->vals[t_id];
      cs_real_t *v_r12 = CS_F_(r12)->vals[t_id];
      cs_real_t *v_r13 = CS_F_(r13)->vals[t_id];
      cs_real_t *v_r23 = CS_F_(r23)->vals[t_id];
      cs_real_t *v_eps = CS_F_(eps)->vals[t_id];

      err_sum += _read_turb_array_1d_compat(r, "k", "k", t_id, v_r11);
      err_sum += _read_turb_array_1d_compat(r, "omega", "omega", t_id, v_eps);

      /* Now transform omega to epsilon */

      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_eps[i] = cs_turb_cmu*v_eps[i]*v_r11[i];

      double d2s3 = 2./3.;

      for (cs_lnum_t i = 0; i < n_cells; i++) {
        double d2s3xk = v_r11[i]*d2s3;
        v_r11[i] = d2s3xk;
        v_r22[i] = d2s3xk;
        v_r33[i] = d2s3xk;
        v_r12[i] = 0.;
        v_r13[i] = 0.;
        v_r23[i] = 0.;
      }

      if (err_sum == 0) {
        read_flag[CS_F_(r11)->id] += t_mask;
        read_flag[CS_F_(r22)->id] += t_mask;
        read_flag[CS_F_(r33)->id] += t_mask;
        read_flag[CS_F_(r12)->id] += t_mask;
        read_flag[CS_F_(r23)->id] += t_mask;
        read_flag[CS_F_(r13)->id] += t_mask;
      }
    }

  } else if (itytur_cur == 5) { /* New computation is in v2f; */

    /* if (itytur_old == 2) k-epsilon to v2f (phi-fbar or BL-v2/k)
     *                      k and epsilon already read, keep default
     *                      initializations for phi and f_bar or alpha,
     *                      so nothing special to do here */

    if (itytur_old == 3) { /* Rij to v2f (phi-fbar or BL-v2/k)
                              epsilon alread read; keep default initializations
                              for phi and f_bar or alpha (v2 in phi is not a
                              true Rij component) */

      cs_real_t *v_k = CS_F_(k)->vals[t_id];

      err_sum += _read_turb_array_1d_compat(r, "r11", "R11", t_id, v_k);

      err_sum += _read_turb_array_1d_compat(r, "r22", "R22", t_id, v_tmp);
      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_k[i] += v_tmp[i];

      err_sum += _read_turb_array_1d_compat(r, "r33", "R33", t_id, v_tmp);
      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_k[i] = 0.5 * (v_k[i] + v_tmp[i]);

      if (err_sum == 0)
        read_flag[CS_F_(k)->id] += t_mask;
    }

    /* if (itytur_old == 5) Switch between v2f variants; k, epsilon,
     *                      and phi already read; keep default initialization
     *                      for f_bar or alpha */

    else if (iturb_old == 60) { /* Switch from k-omega to v2f;
                                   k already read; keep default initialization
                                   for f_bar or alpha */

      cs_real_t *v_eps = CS_F_(eps)->vals[t_id];
      const cs_real_t *v_k = CS_F_(k)->vals[t_id];

      err_sum += _read_turb_array_1d_compat(r, "omega", "omega", t_id, v_eps);

      /* Now transform omega to epsilon */

      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_eps[i] = cs_turb_cmu*v_eps[i]*v_k[i];

      if (err_sum == 0)
        read_flag[CS_F_(eps)->id] += t_mask;
    }

  } else if (iturb_cur == 60) { /* New computation is in k-omega */

    if (itytur_old == 2 || iturb_old == 50) { /* k-epsilon or v2f to k-omega;
                                                 k already read */

      cs_real_t *v_omg = CS_F_(omg)->vals[t_id];
      const cs_real_t *v_k = CS_F_(k)->vals[t_id];

      err_sum += _read_turb_array_1d_compat(r, "epsilon", "eps", t_id, v_omg);

      /* Now transform epsilon to omega */

      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_omg[i] /= (cs_turb_cmu*v_k[i]);

      if (err_sum == 0)
        read_flag[CS_F_(omg)->id] += t_mask;

    } else if (itytur_old == 3) { /* Rij to k-omega */

      cs_real_t *v_k = CS_F_(k)->vals[t_id];
      cs_real_t *v_omg = CS_F_(omg)->vals[t_id];

      err_sum += _read_turb_array_1d_compat(r, "r11", "R11", t_id, v_k);

      err_sum += _read_turb_array_1d_compat(r, "r22", "R22", t_id, v_tmp);
      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_k[i] += v_tmp[i];

      err_sum += _read_turb_array_1d_compat(r, "r33", "R33", t_id, v_tmp);
      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_k[i] = 0.5 * (v_k[i] + v_tmp[i]);

      err_sum += _read_turb_array_1d_compat(r, "epsilon", "eps", t_id, v_omg);

      /* Now transform epsilon to omega */

      for (cs_lnum_t i = 0; i < n_cells; i++)
        v_omg[i] /= (cs_turb_cmu*v_k[i]);

      if (err_sum == 0) {
        read_flag[CS_F_(k)->id] += t_mask;
        read_flag[CS_F_(omg)->id] += t_mask;
      }

    }

  } else if (iturb_cur == 70) { /* current computation with the
                                   Spalart Allmaras (SA) model */

    /* TODO perform the conversion from other models to SA. */

  }

  if (err_sum != 0)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error reading turbulence variables from previous model\n"
         "in restart file \"%s\"."),
       cs_restart_get_name(r));

  BFT_FREE(v_tmp);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read field metadata from checkpoint.
 *
 * Old ids associated to each field are determined for future use.
 * Note that when reading legacy files (Code_Saturne version 3.3 and below),
 * the old id will actually be the old scalar id (-1 for others).
 *
 * \param[in, out]  r              associated restart file pointer
 * \param[out]      old_field_map  name to id map of fields in restart file
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_read_field_info(cs_restart_t           *r,
                           cs_map_name_to_id_t  **old_field_map)
{
  cs_lnum_t sizes[2];
  int retcode;

  /* Initialization */

  const int n_fields = cs_field_n_fields();

  *old_field_map = NULL;

  /* Now read field names, in id order */

  cs_lnum_t *type_buf = NULL;
  char *name_buf = NULL;

  retcode = cs_restart_read_section(r,
                                    "fields:sizes",
                                    CS_MESH_LOCATION_NONE,
                                    2,
                                    CS_TYPE_cs_int_t,
                                    sizes);

  if (retcode == CS_RESTART_SUCCESS) {

    const int n_o_fields = sizes[0];
    cs_lnum_t n_old[] = {0, 0, 0, 0};
    cs_lnum_t n_cur[] = {0, 0, 0, 0};

    /* Now read main metadata */

    BFT_MALLOC(name_buf, sizes[1] + 1, char);
    BFT_MALLOC(type_buf, sizes[0], cs_int_t);

    retcode = cs_restart_read_section(r,
                                      "fields:names",
                                      CS_MESH_LOCATION_NONE,
                                      sizes[1],
                                      CS_TYPE_char,
                                      name_buf);

    if (retcode == CS_RESTART_SUCCESS) {

      cs_map_name_to_id_t *o_map = cs_map_name_to_id_create();

      const char *name = name_buf;

      for (int i = 0, j = 0; j < n_o_fields; i++) {
        if (name_buf[i] == '\0') {
          cs_map_name_to_id(o_map, name);
          name = name_buf + i + 1;
          j++;
        }
      }

      *old_field_map = o_map;

      retcode = cs_restart_read_section(r,
                                        "fields:types",
                                        CS_MESH_LOCATION_NONE,
                                        sizes[0],
                                        CS_TYPE_cs_int_t,
                                        type_buf);

      if (retcode != CS_RESTART_SUCCESS) {
        for (int i = 0; i < n_fields; i++)
          type_buf[i] = 0;
      }
    }

    BFT_FREE(name_buf);

    /* Count variables  */

    for (int f_id = 0; f_id < n_fields; f_id++) {
      const cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        if (f->type & CS_FIELD_USER)
          n_cur[1] += 1;
        else
          n_cur[0] += 1;
      }
      else if (f->type & CS_FIELD_PROPERTY) {
        if (f->type & CS_FIELD_USER)
          n_cur[3] += 1;
        else
          n_cur[2] += 1;
      }
    }

    for (int f_id = 0; f_id < sizes[0]; f_id++) {
      if (type_buf[f_id] & CS_FIELD_VARIABLE) {
        if (type_buf[f_id] & CS_FIELD_USER)
          n_old[1] += 1;
        else
          n_old[0] += 1;
      }
      else if (type_buf[f_id] & CS_FIELD_PROPERTY) {
        if (type_buf[f_id] & CS_FIELD_USER)
          n_old[3] += 1;
        else
          n_old[2] += 1;
      }
    }

    /* Warn in case of change */

    if (   n_cur[0] != n_old[0] || n_cur[1] != n_old[1]
        || n_cur[2] != n_old[2] || n_cur[3] != n_old[3])
      bft_printf
        (_("\n"
           "  Warning: the number of variables or properties has been changed\n"
           "           relative to the restart file.\n\n"
           "  current:  %d model and %d user variables\n"
           "            %d model and %d user properties\n"
           "  previous: %d model and %d user variables\n"
           "            %d model and %d user properties\n\n"
           "  The computation continues, with a partial restart.\n"),
         (int)(n_cur[0]), (int)(n_cur[1]), (int)(n_cur[2]), (int)(n_cur[3]),
         (int)(n_old[0]), (int)(n_old[1]), (int)(n_old[2]), (int)(n_old[3]));

    BFT_FREE(type_buf);
  }

  /* Read legacy metadata */

  else if (cs_restart_is_from_ncfd() == 0)
    _read_legacy_field_info(r);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write field metadata to main checkpoint.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_field_info(cs_restart_t  *r)
{
  cs_lnum_t sizes[2];

  int n_fields = cs_field_n_fields();

  sizes[0] = n_fields;
  sizes[1] = 0;

  /* Build field name buffer */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *f = cs_field_by_id(f_id);
    sizes[1] += strlen(f->name) + 1;
  }

  cs_lnum_t  *type_buf;
  char       *name_buf;

  BFT_MALLOC(type_buf, n_fields, cs_lnum_t);
  BFT_MALLOC(name_buf, sizes[1] + 1, char);

  sizes[1] = 0;

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t *f = cs_field_by_id(f_id);

    size_t l = strlen(f->name) + 1;
    memcpy(name_buf + sizes[1], f->name, l);
    sizes[1] += l;

    type_buf[f_id] = f->type;

  }

  /* Now write data */

  cs_restart_write_section(r,
                           "fields:sizes",
                           CS_MESH_LOCATION_NONE,
                           2,
                           CS_TYPE_cs_int_t,
                           sizes);

  cs_restart_write_section(r,
                           "fields:names",
                           CS_MESH_LOCATION_NONE,
                           sizes[1],
                           CS_TYPE_char,
                           name_buf);

  cs_restart_write_section(r,
                           "fields:types",
                           CS_MESH_LOCATION_NONE,
                           n_fields,
                           CS_TYPE_cs_int_t,
                           type_buf);

  BFT_FREE(name_buf);
  BFT_FREE(type_buf);

  bft_printf(_("  Wrote field names and types to checkpoint"
               " at iteration %d: %s\n"),
             cs_glob_time_step->nt_cur, cs_restart_get_name(r));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read variables from checkpoint.
 *
 * \param[in, out]  r              associated restart file pointer
 * \param[in]       old_field_map  name to id map of fields in restart file
 * \param[in]       t_id_flag      -1: all time values; 0: current values;
 *                                 > 0: previous values
 * \param[in, out]  read_flag      optional flag to track fields read,
 *                                 or NULL; set to sum of 2^time_id for fields
 *                                 read (size: n_fields)
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_read_variables(cs_restart_t               *r,
                          const cs_map_name_to_id_t  *old_field_map,
                          int                         t_id_flag,
                          int                         read_flag[])
{
  const int n_fields = cs_field_n_fields();

  /* Initialization */

  int *_read_flag = read_flag;

  if (_read_flag == NULL) {
    BFT_MALLOC(_read_flag, n_fields, int);
    for (int f_id = 0; f_id < n_fields; f_id++)
      _read_flag[f_id] = 0;
  }

  /* Read metadata for turbulent flux models
     (turbulent flux model useful only for warnings, as the turbulent flux
     field id is sufficient  to determine data to be read) */

  _check_turb_flux_model(r, old_field_map);

  /* Read field data */

  for (int f_id = 0; f_id < n_fields; f_id++) {

    /* Base fields */

    const cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      int t_id_s = 0;
      int t_id_e = f->n_time_vals;
      if (t_id_flag == 0)
        t_id_e = 1;
      else if (t_id_flag > 0)
        t_id_s = 1;

      for (int t_id = t_id_s; t_id < t_id_e; t_id++) {
        int t_mask = (t_id == 0) ? 1 : 2 << (t_id-1);
        if (_read_flag[f->id] & t_mask)
          continue;

        int retval = cs_restart_read_field_vals(r, f_id, t_id);

        if (retval == CS_RESTART_SUCCESS)
          _read_flag[f_id] += t_mask;
      }

    }
  }

  /* Read linked field data;
     note that when turbulent fluxes are variables and not properties,
     they should already have been read by the main loop on variables,
     so use of _read_flag here should avoid writing twice.
  */

  cs_restart_read_linked_fields(r,
                                old_field_map,
                                "turbulent_flux_id",
                                _read_flag);

  cs_restart_read_linked_fields(r,
                                old_field_map,
                                "gwf_sorbed_concentration_id",
                                _read_flag);

  cs_restart_read_linked_fields(r,
                                old_field_map,
                                "gwf_precip_concentration_id",
                                _read_flag);

  /* Read and convert turbulence variables in case of model change */

  const int iturb_cur = cs_glob_turb_model->iturb;
  const int iturb_old = _read_turbulence_model(r);

  if (iturb_cur != iturb_old)
    _read_and_convert_turb_variables(r, iturb_cur, iturb_old, 0, _read_flag);

  /* Warning or error for main unread variables */

  if (t_id_flag < 1) {

    if (cs_glob_field_pointers != NULL) {

      /* standard calculation */
      if (CS_F_(p)) {
        if (!(_read_flag[CS_F_(vel)->id] & 1) || !(_read_flag[CS_F_(p)->id] & 1))
          bft_error
            (__FILE__, __LINE__, 0,
             _("Error reading velocity/pressure values in restart file \"%s\"."),
             cs_restart_get_name(r));
        /* ground water flow calculation */
      } else if (CS_F_(head)) {
        if (   !(_read_flag[CS_F_(vel)->id] & 1)
            || !(_read_flag[CS_F_(head)->id] & 1))
          bft_error
            (__FILE__, __LINE__, 0,
             _("Error reading velocity/hydraulic head values in restart "
               "file \"%s\"."),
             cs_restart_get_name(r));
      }

    }

    int n_missing = 0;

    for (int f_id = 0; f_id < n_fields; f_id++) {
      const cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        if (!(_read_flag[f_id] & 1))
          n_missing++;
      }
    }

    if (n_missing > 0) {
      bft_printf
        (_("\n"
           "  Warning: the following variables could not be found or read\n"
           "           in restart file \"%s\", so default initializations\n"
           "           will be used:\n\n"), cs_restart_get_name(r));
      for (int f_id = 0; f_id < n_fields; f_id++) {
        const cs_field_t *f = cs_field_by_id(f_id);
        if (f->type & CS_FIELD_VARIABLE) {
          if (!(_read_flag[f_id] & 1))
            bft_printf("  %s\n", cs_field_get_label(f));
        }
      }
      bft_printf("\n");
    }

  }

  /* Cleanup */

  if (_read_flag != read_flag)
    BFT_FREE(_read_flag);

  bft_printf(_("  Read variables from restart: %s\n"),
             cs_restart_get_name(r));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write variables to checkpoint.
 *
 * \param[in, out]  r          associated restart file pointer
 * \param[in]       t_id_flag  -1: all time values; 0: current values;
 *                             > 0: previous values
 * \param[in, out]  write_flag  optional flag to track fields written, or NULL;
 *                              set to sum of 2^time_id for fields written
 *                              (size: n_fields)
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_variables(cs_restart_t  *r,
                           int            t_id_flag,
                           int            write_flag[])
{
  int n_fields = cs_field_n_fields();

  /* Initialization */

  int *_write_flag = write_flag;

  if (_write_flag == NULL) {
    BFT_MALLOC(_write_flag, n_fields, int);
    for (int f_id = 0; f_id < n_fields; f_id++)
      _write_flag[f_id] = 0;
  }

  /* Write metadata for turbulent flux models
     (turbulent flux model useful only for warnings, as the turbulent flux
     field id is sufficient  to determine data to be written;
     last test on n_turbt is a minor optimization to possibly avoid a call
     to cs_restart_write_linked_fields, but that call would be no
     more costly than the loop to determine the turbulent flux model
     in the first case...) */

  int n_turbt = 0;

  {
    cs_lnum_t *turbt_buf;

    BFT_MALLOC(turbt_buf, n_fields, cs_lnum_t);

    for (int f_id = 0; f_id < n_fields; f_id++)
      turbt_buf[f_id] = 0;

    int k_sca = cs_field_key_id("scalar_id");
    int k_turbt = cs_field_key_id("turbulent_flux_model");

    for (int f_id = 0; f_id < n_fields; f_id++) {
      const cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        int s_num = cs_field_get_key_int(f, k_sca);
        if (s_num > 0) {
          int f_turbt = cs_field_get_key_int(f, k_turbt);
          if (f_turbt > 0) {
            turbt_buf[f_id] = f_turbt;
            n_turbt++;
          }
        }
      }
    }

    if (n_turbt > 0 && t_id_flag < 1)
      cs_restart_write_section(r,
                               "fields:turbulent_flux_model",
                               CS_MESH_LOCATION_NONE,
                               n_fields,
                               CS_TYPE_cs_int_t,
                               turbt_buf);

    BFT_FREE(turbt_buf);
  }

  /* Write field data */

  for (int f_id = 0; f_id < n_fields; f_id++) {

    /* Base fields */

    const cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      int t_id_s = 0;
      int t_id_e = f->n_time_vals;
      if (t_id_flag == 0)
        t_id_e = 1;
      else if (t_id_flag > 0)
        t_id_s = 1;

      for (int t_id = t_id_s; t_id < t_id_e; t_id++) {
        int t_mask = (t_id == 0) ? 1 : 2 << (t_id-1);
        if (_write_flag[f_id] & t_mask)
          continue;

        cs_restart_write_field_vals(r, f_id, t_id);

        _write_flag[f_id] += t_mask;
      }

    }
  }

  /* Write linked field data;
     note that when turbulent fluxes are variables and not properties,
     they should already have been written by the main loop on variables,
     so use of _write_flag here should avoid writing twice.
  */

  if (n_turbt > 0)
    cs_restart_write_linked_fields(r,
                                   "turbulent_flux_id",
                                   _write_flag);

  cs_restart_write_linked_fields(r,
                                 "gwf_sorbed_concentration_id",
                                 _write_flag);

  cs_restart_write_linked_fields(r,
                                 "gwf_precip_concentration_id",
                                 _write_flag);

  /* Cleanup */

  if (_write_flag != write_flag)
    BFT_FREE(_write_flag);

  bft_printf(_("  Wrote main variables to checkpoint: %s\n"),
             cs_restart_get_name(r));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read fields depending on others from checkpoint.
 *
 * Old ids associate to each field are determined for future use.
 * Note that when reading legacy files (Code_Saturne version 3.3 and below),
 * the old id will actually be the old scalar id (-1 for others).
 *
 * \param[in, out]  r              associated restart file pointer
 * \param[in]       old_field_map  name to id map of fields in restart file
 * \param[in]       key            key for field association
 * \param[in, out]  read_flag      optional flag to track fields read, or NULL;
 *                                 set to sum of 2^time_id for fields read,
 *                                 -1 for fields failed to read (size: n_fields)
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_read_linked_fields(cs_restart_t               *r,
                              const cs_map_name_to_id_t  *old_field_map,
                              const char                 *key,
                              int                         read_flag[])
{
  int retcode;

  /* Initialization */

  int n_required = 0, n_legacy_read = 0;

  const int n_fields = cs_field_n_fields();
  const int n_o_fields = cs_map_name_to_id_size(old_field_map);

  const int key_id = cs_field_key_id_try(key);
  const int key_flag = cs_field_key_flag(key_id);
  const int kr = cs_field_key_id_try("restart_name");

  /* First, check if we need to read anything */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (key_flag != 0) {
      if (key_flag == -1 || !(f->type & key_flag))
        continue;
    }
    if (cs_field_get_key_int(f, key_id) > -1)
      n_required += 1;
  }

  if (n_required < 1)
    return;

  /* Now try to read */

  int *_read_flag = read_flag;

  if (_read_flag == NULL) {
    BFT_MALLOC(_read_flag, n_fields, int);
    for (int f_id = 0; f_id < n_fields; f_id++)
      _read_flag[f_id] = 0;
  }

  cs_lnum_t *old_key_val;
  BFT_MALLOC(old_key_val, n_o_fields, cs_lnum_t);

  char *sec_name;
  BFT_MALLOC(sec_name, strlen("fields:") + strlen(key) + 1, char);
  strcpy(sec_name, "fields:");
  strcat(sec_name, key);

  /* Read metadata */

  retcode = cs_restart_check_section(r,
                                     sec_name,
                                     CS_MESH_LOCATION_NONE,
                                     n_o_fields,
                                     CS_TYPE_cs_int_t);

  /* Try to read in compatibility mode if section not found */

  if (retcode == CS_RESTART_ERR_EXISTS)
    n_legacy_read = _read_linked_fields_legacy(r, key, _read_flag);

  /* Otherwise try to read section anyways to log warning */

  if (n_legacy_read == 0)
    retcode = cs_restart_read_section(r,
                                      sec_name,
                                      CS_MESH_LOCATION_NONE,
                                      n_o_fields,
                                      CS_TYPE_cs_int_t,
                                      old_key_val);

  BFT_FREE(sec_name);

  if (retcode == CS_RESTART_SUCCESS && n_legacy_read == 0) {

    for (int f_id = 0; f_id < n_fields; f_id++) {

      const cs_field_t *f = cs_field_by_id(f_id);

      if (key_flag != 0) {
        if (key_flag == -1 || !(f->type & key_flag))
          continue;
      }

      const int lnk_f_id = cs_field_get_key_int(f, key_id);

      if (lnk_f_id > -1) {

        cs_field_t *f_lnk = cs_field_by_id(lnk_f_id);
        const char *f_lnk_name = NULL;

        if (_read_flag[lnk_f_id] != 0)
          continue;

        /* Check if dependent field has explicit rename */
        if (kr > -1)
          f_lnk_name = cs_field_get_key_str(f, kr);

        /* Otherwise, determine matching name,
           with possible parent field rename */

        if (f_lnk_name == NULL) {
          const char *f_name = NULL;
          if (kr > -1)
            f_name = cs_field_get_key_str(f, kr);
          if (f_name == NULL)
            f_name = f->name;
          int old_f_id = cs_map_name_to_id_try(old_field_map, f_name);
          if (old_f_id > -1) {
            int old_lnk_f_id = old_key_val[old_f_id];
            if (old_lnk_f_id > -1)
              f_lnk_name
                = cs_map_name_to_id_reverse(old_field_map, old_lnk_f_id);
            else
              f_lnk_name = f_lnk->name;
          }
        }

        /* Now we know which field name to read */

        if (f_lnk_name != NULL) {

          _read_flag[lnk_f_id] = -1;

          for (int t_id = 0; t_id < f_lnk->n_time_vals; t_id++) {
            retcode = _read_field_vals(r, f_lnk_name, t_id, f_lnk);
            if (retcode == CS_RESTART_SUCCESS) {
              if (t_id == 0)
                _read_flag[lnk_f_id] = 1;
              else
                _read_flag[lnk_f_id] += (2 << (t_id-1));
            }
            else
              break;
          }

        }
        else if (_read_flag[lnk_f_id] == 0) {
          _read_flag[lnk_f_id] = -1;
          bft_printf(_("  %s: no matching data for field \"%s\"\n"),
                     cs_restart_get_name(r),  f_lnk->name);
          if (f_lnk_name != NULL && f_lnk_name != f_lnk->name)
            bft_printf(_("      (was named \"%s\")\n"),
                       f_lnk_name);
        }

      }

    }
  }

  BFT_FREE(old_key_val);

  if (read_flag != _read_flag)
    BFT_FREE(_read_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write fields depending on others to checkpoint.
 *
 * \param[in, out]  r           associated restart file pointer
 * \param[in]       key         key for field association
 * \param[in, out]  write_flag  optional flag to track fields written, or NULL;
 *                              set to sum of 2^time_id for fields written
 *                              (size: n_fields)
 * \return  number of fields written
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_write_linked_fields(cs_restart_t  *r,
                               const char    *key,
                               int            write_flag[])
{
  int retval = 0;

  /* Initialization */

  const int n_fields = cs_field_n_fields();

  const int key_id = cs_field_key_id_try(key);
  const int key_flag = cs_field_key_flag(key_id);

  int *_write_flag = write_flag;

  if (_write_flag == NULL) {
    BFT_MALLOC(_write_flag, n_fields, int);
    for (int f_id = 0; f_id < n_fields; f_id++)
      _write_flag[f_id] = 0;
  }

  cs_lnum_t *key_val;
  BFT_MALLOC(key_val, n_fields, cs_lnum_t);

  char *sec_name;
  BFT_MALLOC(sec_name, strlen("fields:") + strlen(key) + 1, char);
  strcpy(sec_name, "fields:");
  strcat(sec_name, key);

  /* Write metadata */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    key_val[f_id] = -1;
    const cs_field_t *f = cs_field_by_id(f_id);
    if (key_flag != 0) {
      if (key_flag == -1 || !(f->type & key_flag))
        continue;
    }
    key_val[f_id] = cs_field_get_key_int(f, key_id);
  }

  cs_restart_write_section(r,
                           sec_name,
                           CS_MESH_LOCATION_NONE,
                           n_fields,
                           CS_TYPE_cs_int_t,
                           key_val);

  BFT_FREE(sec_name);

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const int lnk_f_id = key_val[f_id];
    if (lnk_f_id < 0 || _write_flag[lnk_f_id] != 0)
      continue;

    cs_field_t *f_lnk = cs_field_by_id(lnk_f_id);

    _write_flag[lnk_f_id] = -1;

    /* Now write field values */

    for (int t_id = 0; t_id < f_lnk->n_time_vals; t_id++) {

      cs_restart_write_field_vals(r, lnk_f_id, t_id);

      if (t_id == 0)
        _write_flag[lnk_f_id] = 1;
      else
        _write_flag[lnk_f_id] += (2 << (t_id-1));

    }

    retval += 1;

  }

  BFT_FREE(key_val);

  if (_write_flag != write_flag)
    BFT_FREE(_write_flag);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read boundary condition coefficients for all fields from checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_read_bc_coeffs(cs_restart_t  *r)
{
  int c_id, f_id;

  int errcount = 0;
  const int coupled_key_id = cs_field_key_id_try("coupled");
  const int n_fields = cs_field_n_fields();
  char old_name_xx[128] = "", old_name_yy[128] = "", old_name_zz[128] = "";
  char old_name_xy[128] = "", old_name_yz[128] = "", old_name_xz[128] = "";
  const int kr = cs_field_key_id_try("restart_name");

  /* Loop on all fields, to search for those defined on all cells
     and with BC coefficients */

  for (f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    if (   f->location_id == CS_MESH_LOCATION_CELLS
        && f->bc_coeffs != NULL) {

      /* Check for presence of coefficients */

      int coupled = 0;
      int n_loc_vals = 1;

      int32_t coeff_p[] = {0, 0, 0, 0, 0, 0, 0, 0};

      cs_real_t *p[] = {f->bc_coeffs->a,
                        f->bc_coeffs->b,
                        f->bc_coeffs->af,
                        f->bc_coeffs->bf,
                        f->bc_coeffs->ad,
                        f->bc_coeffs->bd,
                        f->bc_coeffs->ac,
                        f->bc_coeffs->bc};

      for (c_id = 0; c_id < 8; c_id++) {
        if (p[c_id] != NULL) {
          coeff_p[c_id] = 1;
          /* avoid double reads/writes in case of aliasing */
          for (int i = 0; i < c_id; i++) {
            if (p[i] == p[c_id])
              coeff_p[c_id] = 0;
          }
        }
      }

      cs_parall_max(8, CS_INT32, coeff_p);

      if (f->dim > 1 && coupled_key_id > -1)
        coupled = cs_field_get_key_int(f, coupled_key_id);

      for (c_id = 0; c_id < 8; c_id++) {
        int retval;
        char *sec_name = NULL;
        cs_real_t *c = p[c_id];
        const char *name = NULL;
        if (kr > -1)
          name = cs_field_get_key_str(f, kr);
        if (name == NULL)
          name = f->name;
        if (coeff_p[c_id] == 0)
          continue;

        if (coupled) {
          if (c_id %2 == 0)
            n_loc_vals = f->dim;
          else
            n_loc_vals = f->dim * f->dim;
        }
        else { /* uncoupled */
          n_loc_vals = f->dim;
        }

        BFT_MALLOC(sec_name,
                   strlen(name) + strlen(_coeff_name[c_id]) + 3,
                   char);
        sprintf(sec_name, "%s::%s", name, _coeff_name[c_id]);

        retval = cs_restart_check_section(r,
                                          sec_name,
                                          f->location_id,
                                          f->dim,
                                          CS_TYPE_cs_real_t);

        if ( f->dim == 6 && retval == CS_RESTART_ERR_EXISTS){
          sprintf(sec_name, "rij::%s", _coeff_name[c_id]);
            snprintf(old_name_xx, 127, "r11::%s", _coeff_name[c_id]);
            snprintf(old_name_yy, 127, "r22::%s", _coeff_name[c_id]);
            snprintf(old_name_zz, 127, "r33::%s", _coeff_name[c_id]);
            snprintf(old_name_xy, 127, "r12::%s", _coeff_name[c_id]);
            snprintf(old_name_yz, 127, "r23::%s", _coeff_name[c_id]);
            snprintf(old_name_xz, 127, "r13::%s", _coeff_name[c_id]);
          if (c_id %2 == 0) {
            retval =  cs_restart_read_real_6_t_compat(r,
                                                      sec_name,
                                                      old_name_xx,
                                                      old_name_yy,
                                                      old_name_zz,
                                                      old_name_xy,
                                                      old_name_yz,
                                                      old_name_xz,
                                                      f->location_id,
                                                      (cs_real_6_t *)(f->val));
          }
          else {
             retval =  cs_restart_read_real_66_t_compat(r,
                                                        sec_name,
                                                        old_name_xx,
                                                        old_name_yy,
                                                        old_name_zz,
                                                        old_name_xy,
                                                        old_name_yz,
                                                        old_name_xz,
                                                        f->location_id,
                                                        (cs_real_66_t *)(f->val));
          }
        }
        else {
          retval = cs_restart_read_section(r,
                                           sec_name,
                                           3, /* location_id */
                                           n_loc_vals,
                                           CS_TYPE_cs_real_t,
                                           c);
        }
        if (retval != CS_RESTART_SUCCESS)
          errcount += 1;

        BFT_FREE(sec_name);

      } /* End of loop in i (coeff type) */

    } /* End for field with BC coeffs */

  } /* End of loop on fields */

  if (errcount > 0) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("\nSome boundary condition coefficients "
                 "could not be read from a restart file;\n"
                 "they will be initialized with default values.\n\n"));
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write boundary condition coefficients for all fields to checkpoint.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_bc_coeffs(cs_restart_t  *r)
{
  int c_id, f_id;

  const int coupled_key_id = cs_field_key_id_try("coupled");
  const int n_fields = cs_field_n_fields();

  /* Loop on all fields, to search for those defined on all cells
     and with BC coefficients */

  for (f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    if (   f->location_id == CS_MESH_LOCATION_CELLS
        && f->bc_coeffs != NULL) {

      /* Check for presence of coefficients */

      int coupled = 0;
      int n_loc_vals = 1;

      int32_t coeff_p[] = {0, 0, 0, 0, 0, 0, 0, 0};
      cs_real_t *p[] = {f->bc_coeffs->a,
                        f->bc_coeffs->b,
                        f->bc_coeffs->af,
                        f->bc_coeffs->bf,
                        f->bc_coeffs->ad,
                        f->bc_coeffs->bd,
                        f->bc_coeffs->ac,
                        f->bc_coeffs->bc};

      for (c_id = 0; c_id < 8; c_id++) {
        if (p[c_id] != NULL) {
          coeff_p[c_id] = 1;
          /* avoid double reads/writes in case of aliasing */
          for (int i = 0; i < c_id; i++) {
            if (p[i] == p[c_id])
              coeff_p[c_id] = 0;
          }
        }
      }

      cs_parall_max(8, CS_INT32, coeff_p);

      if (f->dim > 1 && coupled_key_id > -1)
        coupled = cs_field_get_key_int(f, coupled_key_id);

      for (c_id = 0; c_id < 8; c_id++) {

        char *sec_name = NULL;

        cs_real_t *c = p[c_id];

        if (coeff_p[c_id] == 0)
          continue;

        if (coupled) {
          if (c_id %2 == 0)
            n_loc_vals = f->dim;
          else
            n_loc_vals = f->dim * f->dim;
        }
        else { /* uncoupled */
          n_loc_vals = f->dim;
        }

        BFT_MALLOC(sec_name,
                   strlen(f->name) + strlen(_coeff_name[c_id]) + 3,
                   char);
        sprintf(sec_name, "%s::%s", f->name, _coeff_name[c_id]);

        cs_restart_write_section(r,
                                 sec_name,
                                 3, /* location_id */
                                 n_loc_vals,
                                 CS_TYPE_cs_real_t,
                                 c);

        BFT_FREE(sec_name);

        if (c != p[c_id])
          BFT_FREE(c);

      } /* End of loop in i (coeff type) */

    } /* End for field with BC coeffs */

  } /* End of loop on fields */

  bft_printf(_("  Wrote boundary condition coefficients to checkpoint: %s\n"),
             cs_restart_get_name(r));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read field values from checkpoint.
 *
 * If the values are not found using the default rules based on the
 * field's name, its name itself, or a "restart_rename" keyed string value,
 * an old name may be used for compatibility with older files.
 * For cell-based fields, the old name base is appended automatically with
 * "_ce_phase01", except for scalars, where the name uses a different scheme,
 * based on "scalaire_ce_%04" % s_num;
 *
 * \param[in, out]  r     associated restart file pointer
 * \param[in]       f_id  field id
 * \param[in]       t_id  time id (0 for current, 1 for previous, ...)
 *
 * \return  CS_RESTART_SUCCESS in case of success, CS_RESTART_ERR_... otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_field_vals(cs_restart_t  *r,
                           int            f_id,
                           int            t_id)
{
  cs_field_t  *f = cs_field_by_id(f_id);
  char sec_name[128], ref_sec_name[128];

  const char *r_name = NULL;

  int retcode = CS_RESTART_SUCCESS;

  /* Check for renaming */

  int kr = cs_field_key_id_try("restart_name");
  if (kr > -1)
    r_name = cs_field_get_key_str(f, kr);

  if (r_name == NULL)
    r_name = f->name;

  /* Check for data; data will be read later, so that compatibility
     checks may be done first; we really try reading the data only
     at the end, so if it is not found, a warning will be logged only
     once (and not once per test), and preferentially use the
     base (non-compatibility) name. */

  if (cs_restart_is_from_ncfd()) {
    if (strcmp(r_name, "velocity") == 0)
      snprintf(sec_name, 127, "%s_1", r_name);
    else if (strcmp(r_name, "enthalpy") == 0)
      snprintf(sec_name, 127, "%s_1", r_name);
    else if (strcmp(r_name, "temperature") == 0)
      snprintf(sec_name, 127, "%s_1", r_name);
    else if (strcmp(r_name, "k") == 0)
      snprintf(sec_name, 127, "%s_1", r_name);
    else if (strcmp(r_name, "epsilon") == 0)
      snprintf(sec_name, 127, "%s_1", r_name);
    else if (strcmp(r_name, "rij") == 0)
      snprintf(sec_name, 127, "%s_1", r_name);
    else if (strcmp(r_name, "pressure") == 0)
      snprintf(sec_name, 127, "%s", r_name);
    else
      snprintf(sec_name, 127, "%s::vals::%d", r_name, t_id);

  } else
    snprintf(sec_name, 127, "%s::vals::%d", r_name, t_id);

  sec_name[127] = '\0';
  strncpy(ref_sec_name, sec_name, 128);

  retcode = cs_restart_check_section(r,
                                     sec_name,
                                     f->location_id,
                                     f->dim,
                                     CS_TYPE_cs_real_t);

  /* Otherwise, try reading with basic (restart) name only if requested */

  if (   (retcode == CS_RESTART_ERR_EXISTS || retcode == CS_RESTART_ERR_N_VALS)
      && r_name != f->name) {
    snprintf(sec_name, 127, "%s", r_name);
    sec_name[127] = '\0';
    retcode = cs_restart_check_section(r,
                                       sec_name,
                                       f->location_id,
                                       f->dim,
                                       CS_TYPE_cs_real_t);

  } else if ((retcode == CS_RESTART_ERR_EXISTS || retcode == CS_RESTART_ERR_N_VALS)
             && cs_restart_is_from_ncfd() ) {
    /* If restart from NCFD, ensure that we can read scalar fields! */
    snprintf(sec_name, 127, "%s", r_name);
    sec_name[127] = '\0';
    retcode = cs_restart_check_section(r,
                                       sec_name,
                                       f->location_id,
                                       f->dim,
                                       CS_TYPE_cs_real_t);
  }

  /* Read data if found */

  if ((cs_restart_is_from_ncfd() && strcmp(f->name, "pressure") != 0) ||
      cs_restart_is_from_ncfd() == 0 ) {

    if (retcode == CS_RESTART_SUCCESS)
      retcode = cs_restart_read_section(r,
                                        sec_name,
                                        f->location_id,
                                        f->dim,
                                        CS_TYPE_cs_real_t,
                                        f->vals[t_id]);

    /* Try reading in compatibility mode if not found */

    else if (   retcode == CS_RESTART_ERR_EXISTS
             || retcode == CS_RESTART_ERR_N_VALS) {

      retcode = _read_field_vals_legacy(r, r_name, t_id, f);

      /* if data is still not found, try reading in normal mode,
         so as to have warning/error messages relative to the
         current naming scheme. */

      if (   retcode == CS_RESTART_ERR_EXISTS
          || retcode == CS_RESTART_ERR_N_VALS)
        retcode = cs_restart_read_section(r,
                                          ref_sec_name,
                                          f->location_id,
                                          f->dim,
                                          CS_TYPE_cs_real_t,
                                          f->vals[t_id]);

    }

  } else {
    /* Special case for pressure:
     * In NCFD it is the total pressure which is stored. Hence we need
     * to compute P = P_tot - Phydro - (P0-Pred0)
     */
      retcode = cs_restart_read_section(r,
                                        sec_name,
                                        f->location_id,
                                        f->dim,
                                        CS_TYPE_cs_real_t,
                                        f->vals[t_id]);

      snprintf(sec_name, 127, "%s", "hydro_pressure");
      sec_name[127] = '\0';
      const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
      cs_real_t *v_tmp;
      BFT_MALLOC(v_tmp, n_cells, cs_real_t);

      retcode = cs_restart_read_section(r,
                                        sec_name,
                                        f->location_id,
                                        f->dim,
                                        CS_TYPE_cs_real_t,
                                        v_tmp);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        f->vals[t_id][c_id] -= v_tmp[c_id]
                             + cs_glob_fluid_properties->p0
                             - cs_glob_fluid_properties->pred0;
      }

      BFT_FREE(v_tmp);
  }

  return retcode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write field values to checkpoint.
 *
 * \param[in, out]  r     associated restart file pointer
 * \param[in]       f_id  field id
 * \param[in]       t_id  time id (0 for current, 1 for previous, ...)
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_field_vals(cs_restart_t  *r,
                            int            f_id,
                            int            t_id)
{
  cs_field_t  *f = cs_field_by_id(f_id);
  char sec_name[128];

  snprintf(sec_name, 127, "%s::vals::%d", f->name, t_id);

  cs_restart_write_section(r,
                           sec_name,
                           f->location_id,
                           f->dim,
                           CS_TYPE_cs_real_t,
                           f->vals[t_id]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read restart time step info.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_read_time_step_info(cs_restart_t  *r)
{
  int retval;
  cs_int_t _n_ts = -1;
  cs_real_t _ts = -1;

  /* First syntax */

  retval = cs_restart_read_section(r,
                                   "nbre_pas_de_temps",
                                   0,
                                   1,
                                   CS_TYPE_cs_int_t,
                                   &_n_ts);
  if (retval == CS_RESTART_SUCCESS)
    retval = cs_restart_read_section(r,
                                     "instant_precedent",
                                     0,
                                     1,
                                     CS_TYPE_cs_real_t,
                                     &_ts);

  /* Second syntax */

  else {
    retval = cs_restart_read_section(r,
                                     "ntcabs",
                                     0,
                                     1,
                                     CS_TYPE_cs_int_t,
                                     &_n_ts);
    if (retval == CS_RESTART_SUCCESS)
      retval = cs_restart_read_section(r,
                                       "ttcabs",
                                       0,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       &_ts);
  }

  if (retval == CS_RESTART_SUCCESS)
    cs_time_step_define_prev(_n_ts, _ts);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Loop over all fields and save them in the restart file which id is
 *        passed in argument if it matches their "restart_file" key value.
 *
 * \param[in, out]  r        associated restart file pointer
 * \param[in]       r_id     value of the key "restart_file"
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_fields(cs_restart_t        *r,
                        cs_restart_file_t    r_id)
{
  int n_fields = cs_field_n_fields();
  const int restart_file_key_id = cs_field_key_id("restart_file");

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, restart_file_key_id) == r_id) {
      cs_restart_write_field_vals(r, f_id, 0);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Loop over all fields and read them in the restart file which id is
 *        passed in argument if it matches their "restart_file" key value.
 *
 * \param[in, out]  r        associated restart file pointer
 * \param[in]       r_id     value of the key "restart_file"
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_read_fields(cs_restart_t       *r,
                       cs_restart_file_t   r_id)
{
  int n_fields = cs_field_n_fields();
  const int restart_file_key_id = cs_field_key_id("restart_file");

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, restart_file_key_id) == r_id) {
      cs_restart_read_field_vals(r, f_id, 0);
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
