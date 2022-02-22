/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_DLOPEN)
#include <dlfcn.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal_extract.h"

#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_equation_iterative_solve.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_post.h"
#include "cs_restart.h"
#include "cs_selector.h"
#include "cs_thermal_model.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atmo.h"
#include "cs_atmo_aerosol.h"
#include "cs_atmo_aerosol_ssh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Pointer to SSH-aerosol .so */

static void *_aerosol_so = NULL;
static const char _lib_path[] = "libssh-aerosol.so";

static bool _allow_ssh_postprocess = false;
static bool _update_ssh_thermo = false;
static bool _verbose = false;
static cs_real_t _ssh_time_offset = 0.0;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_DLOPEN)

/*----------------------------------------------------------------------------
 * Call a function of the shared library
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol
 *
 *----------------------------------------------------------------------------*/

static void
_call(void               *handle,
      const char         *name)
{
  typedef void* (*_tmp_sshaerosol_t)(void);
  _tmp_sshaerosol_t fct
    = (_tmp_sshaerosol_t) cs_base_get_dl_function_pointer(handle,
                                                          name,
                                                          true);
  fct();
}

/*----------------------------------------------------------------------------
 * Send a boolean to SSH-aerosol
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in SSH-aerosol
 *   flag             --> boolean exchanged with the external code
 *----------------------------------------------------------------------------*/

static void
_send_bool(void               *handle,
           const char         *name,
           bool                flag)
{

  typedef void* (*_tmp_sshaerosol_t)(bool*);
  _tmp_sshaerosol_t fct =
    (_tmp_sshaerosol_t) cs_base_get_dl_function_pointer(handle,
                                                        name,
                                                        true);
  fct(&flag);
}

/*----------------------------------------------------------------------------
 * Send a double to SSH-aerosol
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in SSH-aerosol
 *   val              --> double exchanged with the external code
 *
 *----------------------------------------------------------------------------*/

static void
_send_double(void               *handle,
             const char         *name,
             cs_real_t           val)
{
  typedef void* (*_tmp_sshaerosol_t)(double*);
  _tmp_sshaerosol_t fct
    = (_tmp_sshaerosol_t) cs_base_get_dl_function_pointer(handle,
                                                          name,
                                                          true);
  double tmp = val;
  fct(&tmp);
}

/*----------------------------------------------------------------------------
 * Receive a boolean from SSH-aerosol, returns it
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in SSH-aerosol
 *
 * return: the boolean received from the external code
 *----------------------------------------------------------------------------*/

static bool
_recv_bool(void               *handle,
           const char         *name)
{
  typedef bool (*_tmp_sshaerosol_t)(void);
  _tmp_sshaerosol_t fct
    = (_tmp_sshaerosol_t) cs_base_get_dl_function_pointer(handle,
                                                          name,
                                                          true);
  bool res = fct();

  return res;
}

/*----------------------------------------------------------------------------
 * Receive a int from SSH-aerosol, returns it
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in SSH-aerosol
 *
 * return: the integer received from the external code
 *----------------------------------------------------------------------------*/

static int
_recv_int(void               *handle,
          const char         *name)
{
  typedef int (*_tmp_sshaerosol_t)(void);
  _tmp_sshaerosol_t fct
    = (_tmp_sshaerosol_t)cs_base_get_dl_function_pointer(handle,
                                                         name,
                                                         true);
  int res = fct();

  return res;
}

/*----------------------------------------------------------------------------
 * Receive a double from SSH-aerosol, return it
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in SSH-aerosol
 *
 * return: the double received from the external code
 *----------------------------------------------------------------------------*/

static cs_real_t
_recv_double(void               *handle,
             const char         *name)
{
  typedef double (*_tmp_sshaerosol_t)(void);
  _tmp_sshaerosol_t fct
    = (_tmp_sshaerosol_t)cs_base_get_dl_function_pointer(handle,
                                                         name,
                                                         true);
  cs_real_t res = fct();

  return res;
}

/*----------------------------------------------------------------------------
 * Send a char array to SSH-aerosol
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in SSH-aerosol
 *   array            <-- array exchanged with the external code
 *----------------------------------------------------------------------------*/

static void
_exchange_char_array(void               *handle,
                     const char         *name,
                     const char         *array)
{
  typedef void* (*_tmp_sshaerosol_t)(const char*);
  _tmp_sshaerosol_t fct
    = (_tmp_sshaerosol_t) cs_base_get_dl_function_pointer(handle,
                                                          name,
                                                          true);
  fct(array);
}

/*----------------------------------------------------------------------------
 * Read the name of given aerosol from SSH-aerosol
 *
 * parameters:
 *   id               <-- id of the aerosol
 *   name             --> name of the aerosol
 *----------------------------------------------------------------------------*/

static void
_sshaerosol_get_aero_name(const int *id, char *name)
{
  typedef void* (*_tmp_sshaerosol_t)(const int*, char *);
  _tmp_sshaerosol_t fct
    = (_tmp_sshaerosol_t) cs_base_get_dl_function_pointer(_aerosol_so,
                                                          "api_sshaerosol_get_aero_name_",
                                                          true);
  fct(id, name);
}

/*----------------------------------------------------------------------------
 * Exchange a double array with SSH-aerosol
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in SSH-aerosol
 *   array            <-> array exchanged with the external code
 *----------------------------------------------------------------------------*/

static void
_exchange_double_array(void               *handle,
                       const char         *name,
                       cs_real_t          *array)
{
  typedef void* (*_tmp_sshaerosol_t)(double*);
  _tmp_sshaerosol_t fct
    = (_tmp_sshaerosol_t)cs_base_get_dl_function_pointer(handle,
                                                         name,
                                                         true);
  fct(array);
}

#endif /* defined(HAVE_DLOPEN)*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function initializes SSH-aerosol.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_ssh_initialize(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH);

  cs_atmo_chemistry_t *at_chem = cs_glob_atmo_chemistry;

#if defined(HAVE_DLOPEN)
  /* Load the shared library */
  if (_verbose)
    bft_printf(" Initialize shared library for aerosol chemistry:\n    %s \n",
               _lib_path);
  _aerosol_so = cs_base_dlopen(_lib_path);

  /* Declare SSH-aerosol as not running standalone */
  _send_bool(_aerosol_so,
             "api_sshaerosol_set_standalone_",
             false);
  if (_verbose) bft_printf(" Set sshaerosol standalone to false.\n");

  /* Force SSH-aerosol to write output to a log file.
   * Only on rank 0 in parallel */
  _send_bool(_aerosol_so,
             "api_sshaerosol_set_logger_",
             (cs_glob_rank_id <= 0) ? true : false);
  if (_verbose) bft_printf(" Set sshaerosol logger to true on rank master.\n");

  /* Initialize SSH-aerosol (default file name is namelist.ssh) */
  {
    const int _namelist_len = 401;
    char namelist_ssh[_namelist_len];
    for (int i = 0; i < _namelist_len; i++) namelist_ssh[i] = '\0';
    if (at_chem->aero_file_name == NULL) {
      strcpy(namelist_ssh, "namelist.ssh");
    } else {
      strcpy(namelist_ssh, at_chem->aero_file_name);
    }
    _exchange_char_array(_aerosol_so,
                         "api_sshaerosol_initialize_",
                         &namelist_ssh[0]);
    _call(_aerosol_so, "api_sshaerosol_init_distributions_");
    if (_verbose) bft_printf(" Shared library sshaerosol initialized.\n");
  }

  /* Using this is not recommended */
  if (_allow_ssh_postprocess && cs_glob_rank_id <= 0) {

    /* InitOutput */
    _call(_aerosol_so, "api_sshaerosol_initoutput_");

    /* Report */
    _call(_aerosol_so, "api_sshaerosol_report_");

    /* Output */
    _call(_aerosol_so, "api_sshaerosol_output_");

  }

  /* If homogeneous time step => set initial_time and set time step */
  if (   cs_glob_time_step_options->idtvar == CS_TIME_STEP_CONSTANT
      || cs_glob_time_step_options->idtvar == CS_TIME_STEP_ADAPTIVE) {

    /* This is used and saved: time in code_saturne starts at zero */
    {
      _ssh_time_offset = _recv_double(_aerosol_so,
                                      "api_sshaerosol_get_initial_t_");
      if (_verbose)
        bft_printf(" Initial time from SSH-aerosol: %f\n", _ssh_time_offset);
    }

    /* Grab initial time and time step from code_saturne */
    /* FIXME: this is not the initial time read in the meteo / chemistry file */
    cs_real_t initial_time = cs_glob_time_step->t_cur + _ssh_time_offset;
    cs_real_t dt = (cs_glob_time_step_options->idtvar == CS_TIME_STEP_ADAPTIVE) ?
                    CS_F_(dt)->val[0] : cs_glob_time_step->dt_ref;

    /* Set initial time, current time and time step in SSH */
    _send_double(_aerosol_so, "api_sshaerosol_set_initial_t_", initial_time);
    _send_double(_aerosol_so, "api_sshaerosol_set_current_t_", initial_time);
    _send_double(_aerosol_so, "api_sshaerosol_set_dt_", dt);

  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Time scheme currently incompatible with SSH-aerosol\n"));
  }

  /* InitPhoto */
  if (at_chem->chemistry_with_photolysis)
    _call(_aerosol_so, "api_sshaerosol_initphoto_");

  /* Last safety check */
  if (  _recv_bool(_aerosol_so, "api_sshaerosol_get_logger_")
      && cs_glob_rank_id > 0)
    bft_printf(" Warning: SSH-logger is not parallel.\n");

  /***
   * At this stage, the external aerosol code SSH-aerosol is fully initialized
   *
   * Below, we initialize some code_saturne specific structures / options
   ***/

  /* Get the number of aerosol layers */
  at_chem->n_layer = _recv_int(_aerosol_so,
                               "api_sshaerosol_get_n_aerosol_layers_");

  /* Get the number of aerosols */
  at_chem->n_size = _recv_int(_aerosol_so,
                              "api_sshaerosol_get_nsize_");

  /* Use shorter names for clarity */
  const int nsp = at_chem->n_species;
  const int nlr = at_chem->n_layer;
  const int nsz = at_chem->n_size;

  /* Reallocate arrays */
  BFT_REALLOC(at_chem->species_to_field_id, nsp + nsz * (nlr + 1), int);
  BFT_REALLOC(at_chem->species_to_scalar_id, nsp + nsz * (nlr + 1), int);

  /* For all aerosols */
  for (int i = nsp; i < nsp + nsz * (nlr + 1);  i++ ) {

    /* Build the name of the field */
    char name[512] = "";

    /* Number of the layer [1, N_layer + 1] */
    int ilr = 1 + (i - nsp) / nsz;

    /* Number of the aerosol [1, N_size] */
    const int isize = 1 + (i - nsp) - (ilr - 1) * nsz;

    /* Get the prefix */
    if (ilr <= at_chem->n_layer) {
      /* If possible, import the name from SSH */
      if (1 == _recv_int(_aerosol_so,
                         "api_sshaerosol_get_nlayer_")) {
        char ssh_name[81];
        _sshaerosol_get_aero_name(&ilr, &ssh_name[0]);
        snprintf(name, 81, "%s", ssh_name);
      }
      else {
        if (ilr < 0)
          bft_error(__FILE__, __LINE__, 0,
                    _("Atmospheric aerosols: Number of layers negative."));
        if (ilr > 9999)
          bft_error(__FILE__, __LINE__, 0,
                    _("Atmospheric aerosols: Number of layers above limit."));
        sprintf(name, "aerosol_layer_%04d", ilr);
      }
    }
    else {
      strcpy(name, "aerosol_num");
    }

    /* Get the suffix */
    if (isize < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Atmospheric aerosols : Number of aerosols negative."));
    if (isize > 999)
      bft_error(__FILE__, __LINE__, 0,
                _("Atmospheric aerosols : Number of aerosols above limit."));
    char suffix[5];
    sprintf(suffix, "_%03d", isize);
    strcat(name, suffix);

    /* Field of dimension 1 */
    at_chem->species_to_field_id[i]
      = cs_variable_field_create(name, name, CS_MESH_LOCATION_CELLS, 1);

    /* Scalar field, store in isca_chem/species_to_scalar_id (FORTRAN/C) array */
    at_chem->species_to_scalar_id[i]
      = cs_add_model_field_indexes(at_chem->species_to_field_id[i]);

  }

#else
  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to initialize: %s\n"), _lib_path);
#endif

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function finalizes SSH-aerosol.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_ssh_finalize(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH);

#if defined(HAVE_DLOPEN)
  /* Finalize SSH */
  _call(_aerosol_so, "api_sshaerosol_finalize_");

  /* dlclose: release the shared library */
  cs_base_dlclose(_lib_path, _aerosol_so);
#else
  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to close: %s\n"), _lib_path);
#endif

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function uses the given array to update the aerosol
 *        concentrations and numbers in SSH-aerosol.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_ssh_set_aero(cs_real_t* array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH);

#if defined(HAVE_DLOPEN)
  const int _size = cs_glob_atmo_chemistry->n_layer
                  * cs_glob_atmo_chemistry->n_size;

  /* Set the aerosols concentrations */
  _exchange_double_array(_aerosol_so,
                         "api_sshaerosol_set_aero_",
                         array);

  /* Set the aerosols numbers */
  _exchange_double_array(_aerosol_so,
                         "api_sshaerosol_set_aero_num_",
                         &(array[_size]));
#else
  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to use %s\n"), _lib_path);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with aerosol concentrations and
 *        numbers from SSH-aerosol.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_ssh_get_aero(cs_real_t* array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH);

#if defined(HAVE_DLOPEN)
  const int _size = cs_glob_atmo_chemistry->n_layer
                  * cs_glob_atmo_chemistry->n_size;

  /* Get the aerosols concentrations */
  {
    double data[_size];
    _exchange_double_array(_aerosol_so,
                           "api_sshaerosol_get_aero_",
                           (double *)&data);
    for (int i = 0; i < _size; i++)
      array[i] = data[i];
  }

  /* Get the aerosols numbers */
  {
    double data[cs_glob_atmo_chemistry->n_size];
    _exchange_double_array(_aerosol_so,
                           "api_sshaerosol_get_aero_num_",
                           (double *)&data);
    for (int i = 0; i < cs_glob_atmo_chemistry->n_size; i++)
      array[_size + i] = data[i];
  }
#else
  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to use %s\n"), _lib_path);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function uses the given array to update the gas concentrations
 *        in SSH-aerosol.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_ssh_set_gas(cs_real_t* array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH);

#if defined(HAVE_DLOPEN)
  /* Set the gas concentrations */
  _exchange_double_array(_aerosol_so,
                         "api_sshaerosol_set_gas_",
                         array);
#else
  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to use %s\n"), _lib_path);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with gas concentrations from
 *        SSH-aerosol.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_ssh_get_gas(cs_real_t* array)
{
  assert(cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH);

#if defined(HAVE_DLOPEN)
  /* Get the gas concentrations */
  double data[cs_glob_atmo_chemistry->n_species];
  _exchange_double_array(_aerosol_so,
                         "api_sshaerosol_get_gas_",
                         (double *)&data);
  for (int i = 0; i < cs_glob_atmo_chemistry->n_species; i++)
    array[i] = data[i];
#else
  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to use %s\n"), _lib_path);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes a time step of gaseous chemistry and aerosols
 *        dynamic using SSH-aerosol.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_ssh_time_advance(void)
{
  assert(cs_glob_atmo_chemistry->aerosol_model == CS_ATMO_AEROSOL_SSH);

#if defined(HAVE_DLOPEN)
  const cs_domain_t *domain = cs_glob_domain;
  const cs_mesh_t *m = domain->mesh;

  /* If homogeneous time step => set current_time and set time step
     and update photolysis*/
  if (cs_glob_time_step_options->idtvar == CS_TIME_STEP_CONSTANT
      || cs_glob_time_step_options->idtvar == CS_TIME_STEP_ADAPTIVE) {

    double dt = (cs_glob_time_step_options->idtvar == CS_TIME_STEP_ADAPTIVE) ?
                 CS_F_(dt)->val[0] : cs_glob_time_step->dt_ref;
    double current_time = cs_glob_time_step->t_cur + _ssh_time_offset - dt;

    /* Set the current time */
    _send_double(_aerosol_so, "api_sshaerosol_set_current_t_", current_time);

    /* Set the time step */
    _send_double(_aerosol_so, "api_sshaerosol_set_dt_", dt);

    /* Update the photolysis if needed */
    if (cs_glob_atmo_chemistry->chemistry_with_photolysis)
      _call(_aerosol_so, "api_sshaerosol_updatephoto_");

  } else {
    bft_error(__FILE__, __LINE__, 0,
              _("Time scheme currently incompatible with SSH-aerosol\n"));
  }

  /* Loop over cells, update chemistry and aerosols */
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++)
  {

    /* Conversion from ppm to microg / m^3 */
    const cs_real_t ppm_to_microg = 1e-3 * CS_F_(rho)->val[cell_id];
    const cs_real_t microg_to_ppm = 1. / ppm_to_microg;

    /* Set the Pressure */
    if (_update_ssh_thermo) {
      double pres = cs_field_by_name("total_pressure")->val[cell_id];
      _send_double(_aerosol_so, "api_sshaerosol_set_pressure_", pres);
    }

    /* Set the Temperature (K) */
    if (_update_ssh_thermo) {
      double temp;
      if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE) {
        temp = cs_field_by_name("temperature")->val[cell_id];
        if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS)
          temp -= cs_physical_constants_celsius_to_kelvin;
      } else {
        /* We have enthalpy, use reference temperature */
        temp = cs_glob_fluid_properties->t0;
      }
      _send_double(_aerosol_so, "api_sshaerosol_set_temperature_", temp);
    }

    /* Set the pH */
    if (false) {
      /* Atmospheric module of code_saturne does not provide pH yet */
      double ph = 7;
      _send_double(_aerosol_so, "api_sshaerosol_set_ph_", ph);
    }

    /* Set the relative humidity */
    if (_update_ssh_thermo) {
      cs_field_t* fld = cs_field_by_name_try("total_water");
      if (fld != NULL) {
        cs_real_t totwt = fld->val[cell_id];
        cs_real_t liqwt = cs_field_by_name("liquid_water")->val[cell_id];
        if (fabs(1. - liqwt) < cs_math_epzero)
          bft_error
            (__FILE__,__LINE__, 0,
             _("Error when computing the relative humidity for SSH-aerosol."));
        double rh = (totwt - liqwt)/(1. - liqwt);
        _send_double(_aerosol_so, "api_sshaerosol_set_relhumidity_", rh);
      }
    }

    /* Update the specific humidity */
    if (_update_ssh_thermo)
      _call(_aerosol_so, "api_sshaerosol_update_humidity_");

    /* Update SSH gaseous concentrations */
    {
      double data[cs_glob_atmo_chemistry->n_species];
      for (int i = 0; i < cs_glob_atmo_chemistry->n_species; i++) {
        const int fid = cs_glob_atmo_chemistry->species_to_field_id[i];
        data[i] = cs_field_by_id(fid)->val[cell_id] * ppm_to_microg;
      }
      cs_atmo_aerosol_ssh_set_gas(data);
    }

    /* Update SSH aerosols concentrations and numbers */
    {
      const int _size = cs_glob_atmo_chemistry->n_layer
                      * cs_glob_atmo_chemistry->n_size;
      const int _sizetot = _size + cs_glob_atmo_chemistry->n_size;
      double data[_sizetot];

      /* Concentrations are converted to microg / m^3 */
      for (int i = 0; i < _size; i++) {
        const int ics = i + cs_glob_atmo_chemistry->n_species;
        const int fid = cs_glob_atmo_chemistry->species_to_field_id[ics];
        data[i] = cs_field_by_id(fid)->val[cell_id] * ppm_to_microg;
      }
      /* Numbers are converted to molecules / m^3 */
      for (int i = 0; i < cs_glob_atmo_chemistry->n_size; i++) {
        const int ii = i + _size;
        const int ics = ii + cs_glob_atmo_chemistry->n_species;
        const int fid = cs_glob_atmo_chemistry->species_to_field_id[ics];
        data[ii] = cs_field_by_id(fid)->val[cell_id] * CS_F_(rho)->val[cell_id];
      }
      cs_atmo_aerosol_ssh_set_aero(data);
    }

    /* Update concentration-dependent arrays in SSH-aerosol */
    _call(_aerosol_so, "api_sshaerosol_init_again_");

    /* Emissions */
    _call(_aerosol_so, "api_sshaerosol_emission_");

    /* Call the gaseous chemistry */
    _call(_aerosol_so, "api_sshaerosol_gaschemistry_");

    /* Call the aerosols dynamic */
    _call(_aerosol_so, "api_sshaerosol_aerodyn_");

    /* Using this is not recommended */
    if (_allow_ssh_postprocess && cs_glob_rank_id <= 0 && cell_id == 0) {
      _call(_aerosol_so, "api_sshaerosol_output_");
    }

    /* Update CS gaseous concentrations */
    if (!cs_glob_atmo_chemistry->frozen_gas_chem) {
      double data[cs_glob_atmo_chemistry->n_species];
      cs_atmo_aerosol_ssh_get_gas(data);
      for (int i = 0; i < cs_glob_atmo_chemistry->n_species; i++) {
        const int fid = cs_glob_atmo_chemistry->species_to_field_id[i];
        cs_field_by_id(fid)->val[cell_id] = data[i] * microg_to_ppm;
      }
    }

    /* Update the CS aerosols concentrations and numbers */
    {
      const int _size = cs_glob_atmo_chemistry->n_layer
                      * cs_glob_atmo_chemistry->n_size;
      const int _sizetot = _size + cs_glob_atmo_chemistry->n_size;
      double data[_sizetot];

      /* Get concentrations and numbers */
      cs_atmo_aerosol_ssh_get_aero(data);

      /* Update CS aerosols concentrations */
      for (int i = 0; i < _size; i++) {
        const int ics = i + cs_glob_atmo_chemistry->n_species;
        const int fid = cs_glob_atmo_chemistry->species_to_field_id[ics];
        cs_field_by_id(fid)->val[cell_id] = data[i] * microg_to_ppm;
      }

      /* Update CS aerosols numbers */
      for (int i = 0; i < cs_glob_atmo_chemistry->n_size; i++) {
        const int ii = i + _size;
        const int ics = ii + cs_glob_atmo_chemistry->n_species;
        const int fid = cs_glob_atmo_chemistry->species_to_field_id[ics];
        cs_field_by_id(fid)->val[cell_id] = data[ii] / CS_F_(rho)->val[cell_id];
      }

    }

  }

#else
  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to use %s\n"), _lib_path);
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
