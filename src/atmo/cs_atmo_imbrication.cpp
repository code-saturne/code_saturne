/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_measures_util.h"
#include "base/cs_physical_constants.h"
#include "bft/bft_printf.h"

#include "pprt/cs_physical_model.h"
#include "atmo/cs_air_props.h"
#include "atmo/cs_intprf.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo.h"
#include "atmo/cs_atmo_imbrication.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Parameter for "meteo" files
 *----------------------------------------------------------------------------*/

#define line_len 132
#define skip_chars "/#!"

static char **imbrication_files = nullptr;
static const char *imbrication_files_list = nullptr;

static int number_of_files = 0;
// Time sections per files
static int sections_per_file = -1;
// Profile dimension variable
static int thermal_profile_dim = -1;
static int dynamic_profile_dim = -1;

/*----------------------------------------------------------------------------
 * Read data from "meteo" files
 *----------------------------------------------------------------------------*/

// Time variables
static int *years = nullptr;
static int *hours = nullptr;
static int *minutes = nullptr;
static int *ordinals = nullptr;
static cs_real_t *seconds = nullptr;

// Positions
static cs_real_t *xpos = nullptr;
static cs_real_t *ypos = nullptr;
static cs_real_t *ground_pressure = nullptr;
// Vertical grid for temperature and humidity variables
static cs_real_t *zt = nullptr;
static cs_real_t *qw = nullptr;
static cs_real_t *Nc = nullptr;
static cs_real_t *temp_c = nullptr;
// Vertical grid for wind variables
static cs_real_t *zd = nullptr;
static cs_real_t *tke = nullptr;
static cs_real_t *eps = nullptr;
static cs_real_t *vel_u = nullptr;
static cs_real_t *vel_v = nullptr;

/*----------------------------------------------------------------------------
 * Derived data
 *----------------------------------------------------------------------------*/

static cs_real_t *times = nullptr;
static cs_real_t *theta = nullptr;
static cs_real_t *density = nullptr;
static cs_real_t *pressure = nullptr;

/*----------------------------------------------------------------------------
 * Time interpolated profiles
 *----------------------------------------------------------------------------*/

static cs_real_t *ti_vel_u = nullptr;
static cs_real_t *ti_vel_v = nullptr;
static cs_real_t *ti_zt = nullptr;
static cs_real_t *ti_qw = nullptr;
static cs_real_t *ti_Nc = nullptr;
static cs_real_t *ti_zd = nullptr;
static cs_real_t *ti_tke = nullptr;
static cs_real_t *ti_eps = nullptr;
static cs_real_t *ti_temp_c = nullptr;
static cs_real_t *ti_theta = nullptr;
static cs_real_t *ti_density = nullptr;
static cs_real_t *ti_pressure = nullptr;

/*----------------------------------------------------------------------------
 * Additional variables
 *----------------------------------------------------------------------------*/

static cs_real_t *coordinates_th = nullptr;
static cs_real_t *coordinates_dyn = nullptr;
static cs_real_t *influence_param_th = nullptr;
static cs_real_t *influence_param_dyn = nullptr;

/* atmo imbrication options structure */
static cs_atmo_imbrication_t _atmo_imbrication = {
  .imbrication_flag = false,
  .imbrication_verbose = false,
  .cressman_u = false,
  .cressman_v = false,
  .cressman_qw = false,
  .cressman_nc = false,
  .cressman_tke = false,
  .cressman_eps = false,
  .cressman_theta = false,
  .vertical_influence_radius = 100.0,
  .horizontal_influence_radius = 8500.0,
  .id_u     = -1,
  .id_v     = -1,
  .id_qw    = -1,
  .id_nc    = -1,
  .id_tke   = -1,
  .id_eps   = -1,
  .id_theta = -1
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_atmo_imbrication_t *cs_glob_atmo_imbrication = &_atmo_imbrication;

/*============================================================================
 * Private function definitions
 *============================================================================*/

static void
_allocate_all(void)
{
  const int size = sections_per_file*number_of_files;
  const int thermal_size
    = thermal_profile_dim*sections_per_file*number_of_files;
  const int dyn_size = dynamic_profile_dim*sections_per_file*number_of_files;

  cs_atmo_model_t atmo_model_flag
    = static_cast<cs_atmo_model_t>(cs_glob_physical_model_flag[CS_ATMOSPHERIC]);

  CS_MALLOC(hours, size, int);
  CS_MALLOC(years, size, int);
  CS_MALLOC(minutes, size, int);
  CS_MALLOC(ordinals, size, int);
  CS_MALLOC(seconds, size, cs_real_t);

  // --- section of explicit allocations ---
  CS_MALLOC(xpos, size, cs_real_t);
  CS_MALLOC(ypos, size, cs_real_t);
  CS_MALLOC(zt, thermal_size, cs_real_t);
  CS_MALLOC(ground_pressure, size, cs_real_t);

  /* ---------------------------------------------------------------
   * temp_c and qw are always allocated for dry and humid atmosphere
   * --------------------------------------------------------------- */

  if (atmo_model_flag != CS_ATMO_OFF) {
    CS_MALLOC(qw, thermal_size, cs_real_t);
    CS_MALLOC(temp_c, thermal_size, cs_real_t);
  }

  /* --------------------------------------------------------------
   * Nc only used in humid atmosphere
   * -------------------------------------------------------------- */

  if (atmo_model_flag == CS_ATMO_HUMID)
    CS_MALLOC(Nc, thermal_size, cs_real_t);

  CS_MALLOC(vel_u, dyn_size, cs_real_t);
  CS_MALLOC(vel_v, dyn_size, cs_real_t);
  CS_MALLOC(zd, dyn_size, cs_real_t);
  CS_MALLOC(tke, dyn_size, cs_real_t);
  CS_MALLOC(eps, dyn_size, cs_real_t);
}

/*----------------------------------------------------------------------------
 * \brief Search for the position of a value in an array,
 *        assuming that the array is sorted in a strictly increasing order
 *
 * return if possible lower,upper such that :
 *                    the_array(lower) <= the_value <= the_array(upper)
 *        otherwise :
 *          lower==upper if the_value<the_array(first)
 *                       or if the_value>the_array(last)
 *          lower> upper if none of the previous cases applies (anomaly)
 *
 * \param[in]   the_array   an array
 * \param[in]   the_value   a particular value
 * \param[out]  lower       the index of the first membre of the array
 *                          lower than the value
 * \param[out]   upper      the index of the first membre of the array
 *                          greater than the value
 *----------------------------------------------------------------------------*/

static void
_get_index(const cs_real_t the_array[],
           cs_real_t the_value,
           int &lower,
           int &upper)
{
  int dmin = 0;
  int dmax = sections_per_file - 1;
  for (int i = dmin; i < dmax; i++) {
    if (the_array[i] <= the_value && the_value <= the_array[i+1]) {
      lower = i;
      upper = i+1;
      return;
    }
  }

  if (the_value < the_array[dmin]) {
    lower = dmin;
    upper = dmin;
    return;
  }

  if (the_value > the_array[dmax]) {
    lower = dmax;
    upper = dmax;
    return;
  }

  lower = dmax;
  upper = dmin;
}

/*----------------------------------------------------------------------------
 * \brief Interpolates a "profile" at a given time.
 *        Given a series of profiles varying in time
 *        you get the profile interpolated from them at the given time.
 *
 * \param[in]    the_time                current time
 * \param[in]    the_times               times array
 * \param[in]    the_profiles            input profiles
 * \param[out]   interpolated_profile    output profile
 *----------------------------------------------------------------------------*/

static void
_time_interpolation(const int n_profile,
                    const cs_real_t  the_time,
                    const cs_real_t  the_times[],
                    const cs_real_t  the_profiles[],
                    cs_real_t  interpolated_profile[])
{
  int lower, upper;
  _get_index(the_times, the_time, lower, upper);
  if (lower < upper) {
    const cs_real_t weight = (the_time - the_times[lower]) /
      (the_times[upper] - the_times[lower]);
    if (_atmo_imbrication.imbrication_verbose)
      bft_printf("time_interpolation:: weight=%lf\n", weight);

    for (int i = 0; i < n_profile; i++)
      interpolated_profile[i] =
        the_profiles[i + lower * n_profile] * (1.0 - weight) +
        the_profiles[i + upper * n_profile] * weight;
  }
  else if (lower == upper) {
    for (int i = 0; i < n_profile; i++)
      interpolated_profile[i] = the_profiles[i + lower * n_profile];
  }
  else {
    bft_printf("Time_interpolation: the times array is not increasing\n");
    for (int i = 0; i < sections_per_file; i++)
      bft_printf("Time_interpolation: the_times[%d]=%f\n", i, the_times[i]);
    bft_error(__FILE__, __LINE__, 0,
              _("Time_interpolation stops the calculations."));
  }
}

/*----------------------------------------------------------------------------
 * \brief Time interpolation of all profiles
 *
 * \param[in]   the_time    current time
 *----------------------------------------------------------------------------*/

static void
_interpolate_all_profiles(cs_real_t the_time)
{
  static bool first_call = true;
  static cs_real_t time_was = -1.0;

  if (first_call) {
    time_was = the_time;
    first_call = false;
  }
  else {
    if (fabs(time_was - the_time) < cs_math_epzero) {
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("interpolate_all_profiles:the_time==time_was==%10.05lf\n",
                   the_time);
      return;
    }
  }

  if (_atmo_imbrication.imbrication_verbose)
    bft_printf("interpolate_all_profiles:the_time==%10.05lf\n"
               "interpolate_all_profiles:time_was==%10.05lf\n"
               , the_time, time_was);

  time_was = the_time;

  /* zt */
  static int ti_zt_size = 0;
  if (zt != nullptr) {
    int lb = 0;
    int ub = thermal_profile_dim - 1;
    int required_size = (ub - lb + 1) * number_of_files;

    if (ti_zt_size != required_size && ti_zt != nullptr) {
      CS_FREE(ti_zt);
      ti_zt = nullptr;
      ti_zt_size = 0;
    }

    if (ti_zt == nullptr) {
      CS_MALLOC(ti_zt, required_size, cs_real_t);
      ti_zt_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_zt = ti_zt + i * (ub - lb + 1);
      cs_real_t *_zt = zt + i * thermal_profile_dim * sections_per_file;
      _time_interpolation(thermal_profile_dim,
                          the_time,
                          _times,
                          _zt,
                          _ti_zt);
    }
  }

  /* temp_c */
  static int ti_temp_c_size = 0;
  if (temp_c != nullptr) {
    int lb = 0;
    int ub = thermal_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_temp_c_size != required_size && ti_temp_c != nullptr) {
      CS_FREE(ti_temp_c);
      ti_temp_c = nullptr;
      ti_temp_c_size = 0;
    }

    if (ti_temp_c == nullptr) {
      CS_MALLOC(ti_temp_c, required_size, cs_real_t);
      ti_temp_c_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_temp_c = ti_temp_c + i * dim1;
      cs_real_t *_temp_c = temp_c + i * thermal_profile_dim * sections_per_file;
      _time_interpolation(dim1,
                          the_time,
                          _times,
                          _temp_c,
                          _ti_temp_c);
    }
  }

  /* qw */
  static int ti_qw_size = 0;
  if (qw != nullptr) {
    int lb = 0;
    int ub = thermal_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_qw_size != required_size && ti_qw != nullptr) {
      CS_FREE(ti_qw);
      ti_qw = nullptr;
      ti_qw_size = 0;
    }

    if (ti_qw == nullptr) {
      CS_MALLOC(ti_qw, required_size, cs_real_t);
      ti_qw_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_qw = ti_qw + i * dim1;
      cs_real_t *_qw = qw + i * thermal_profile_dim * sections_per_file;
      _time_interpolation(dim1,
                          the_time,
                          _times,
                          _qw,
                          _ti_qw);
    }
  }

  /* Nc */
  static int ti_Nc_size = 0;
  if (Nc != nullptr) {
    int lb = 0;
    int ub = thermal_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_Nc_size != required_size && ti_Nc != nullptr) {
      CS_FREE(ti_Nc);
      ti_Nc = nullptr;
      ti_Nc_size = 0;
    }

    if (ti_Nc == nullptr) {
      CS_MALLOC(ti_Nc, required_size, cs_real_t);
      ti_Nc_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_Nc = ti_Nc + i * dim1;
      cs_real_t *_Nc = Nc + i * thermal_profile_dim * sections_per_file;
      _time_interpolation(dim1,
                          the_time,
                          _times,
                          _Nc,
                          _ti_Nc);
    }
  }

  /* zd */
  static int ti_zd_size = 0;
  if (zd != nullptr) {
    int lb = 0;
    int ub = dynamic_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_zd_size != required_size && ti_zd != nullptr) {
      CS_FREE(ti_zd);
      ti_zd = nullptr;
      ti_zd_size = 0;
    }

    if (ti_zd == nullptr) {
      CS_MALLOC(ti_zd, required_size, cs_real_t);
      ti_zd_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_zd = ti_zd + i * dim1;
      cs_real_t *_zd = zd + i * dynamic_profile_dim * sections_per_file;
      _time_interpolation(dim1,
                          the_time,
                          _times,
                          _zd,
                          _ti_zd);
    }
  }

  /* vel_u */
  static int ti_vel_u_size = 0;
  if (vel_u != nullptr) {
    int lb = 0;
    int ub = dynamic_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_vel_u_size != required_size && ti_vel_u != nullptr) {
      CS_FREE(ti_vel_u);
      ti_vel_u = nullptr;
      ti_vel_u_size = 0;
    }

    if (ti_vel_u == nullptr) {
      CS_MALLOC(ti_vel_u, required_size, cs_real_t);
      ti_vel_u_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_u = ti_vel_u + i * dim1;
      cs_real_t *_u = vel_u + i * dynamic_profile_dim * sections_per_file;
      _time_interpolation(dim1,
                          the_time,
                          _times,
                          _u,
                          _ti_u);
    }
  }

  /* vel_v */
  static int ti_vel_v_size = 0;
  if (vel_v != nullptr) {
    int lb = 0;
    int ub = dynamic_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_vel_v_size != required_size && ti_vel_v != nullptr) {
      CS_FREE(ti_vel_v);
      ti_vel_v = nullptr;
      ti_vel_v_size = 0;
    }

    if (ti_vel_v == nullptr) {
      CS_MALLOC(ti_vel_v, required_size, cs_real_t);
      ti_vel_v_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_v = ti_vel_v + i * dim1;
      cs_real_t *_v = vel_v + i * dynamic_profile_dim * sections_per_file;
      _time_interpolation(dim1,
                          the_time,
                          _times,
                          _v,
                          _ti_v);
    }
  }

  /* tke */
  static int ti_tke_size = 0;
  if (tke != nullptr) {
    int lb = 0;
    int ub = dynamic_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_tke_size != required_size && ti_tke != nullptr) {
      CS_FREE(ti_tke);
      ti_tke = nullptr;
      ti_tke_size = 0;
    }
    if (ti_tke == nullptr) {
      CS_MALLOC(ti_tke, required_size, cs_real_t);
      ti_tke_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_tke = ti_tke + i * dim1;
      cs_real_t *_tke = tke + i * dynamic_profile_dim * sections_per_file;
      _time_interpolation(dim1, the_time, _times, _tke, _ti_tke);
    }
  }

  /* eps */
  static int ti_eps_size = 0;
  if (eps != nullptr) {
    int lb = 0;
    int ub = dynamic_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_eps_size != required_size && ti_eps != nullptr) {
      CS_FREE(ti_eps);
      ti_eps = nullptr;
      ti_eps_size = 0;
    }
    if (ti_eps == nullptr) {
      CS_MALLOC(ti_eps, required_size, cs_real_t);
      ti_eps_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_eps = ti_eps + i * dim1;
      cs_real_t *_eps = eps + i * dynamic_profile_dim * sections_per_file;
      _time_interpolation(dim1, the_time, _times, _eps, _ti_eps);
    }
  }

  /* pressure */
  static int ti_pressure_size = 0;
  if (pressure != nullptr) {
    int lb = 0;
    int ub = thermal_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_pressure_size != required_size && ti_pressure != nullptr) {
      CS_FREE(ti_pressure);
      ti_pressure = nullptr;
      ti_pressure_size = 0;
    }
    if (ti_pressure == nullptr) {
      CS_MALLOC(ti_pressure, required_size, cs_real_t);
      ti_pressure_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_pressure = ti_pressure + i * dim1;
      cs_real_t *_pressure = pressure + i * thermal_profile_dim * sections_per_file;
      _time_interpolation(dim1, the_time, _times, _pressure, _ti_pressure);
    }
  }

  /* theta */
  static int ti_theta_size = 0;
  if (theta != nullptr) {
    int lb = 0;
    int ub = thermal_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_theta_size != required_size && ti_theta != nullptr) {
      CS_FREE(ti_theta);
      ti_theta = nullptr;
      ti_theta_size = 0;
    }
    if (ti_theta == nullptr) {
      CS_MALLOC(ti_theta, required_size, cs_real_t);
      ti_theta_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_theta = ti_theta + i * dim1;
      cs_real_t *_theta = theta + i * thermal_profile_dim * sections_per_file;
      _time_interpolation(dim1, the_time, _times, _theta, _ti_theta);
    }
  }

  /* density */
  static int ti_density_size = 0;
  if (density != nullptr) {
    int lb = 0;
    int ub = thermal_profile_dim - 1;
    int dim1 = ub - lb + 1;
    int required_size = dim1 * number_of_files;

    if (ti_density_size != required_size && ti_density != nullptr) {
      CS_FREE(ti_density);
      ti_density = nullptr;
      ti_density_size = 0;
    }
    if (ti_density == nullptr) {
      CS_MALLOC(ti_density, required_size, cs_real_t);
      ti_density_size = required_size;
    }

    cs_real_t *_times = times;
    for (int i = 0; i < number_of_files; i++) {
      cs_real_t *_ti_density = ti_density + i * dim1;
      cs_real_t *_density = density + i * thermal_profile_dim * sections_per_file;
      _time_interpolation(dim1, the_time, _times, _density, _ti_density);
    }
  }
}

/*----------------------------------------------------------------------------
 * \brief  Converts a (year,ordinal) date to julian calendar date
 *         for calculating time shifts
 *
 * \param[in]   year
 * \param[in]   ordinal     number of the day in the year
 *              e.g 1st january has ordinal 1
 *                               31 december 365 or 366
 *----------------------------------------------------------------------------*/

static int
_yo2j(const cs_real_t year,
     const cs_real_t ordinal)
{
  int result = ordinal + ( (1461 * (year + 4800)) / 4
                       -    30 - (3 * ((year + 4900) / 100)) / 4
                       +     1 - 32075 ) - 1;

  return result;
}

/*----------------------------------------------------------------------------
 * \brief  Identification of the first and last non white character
 *      of a string
 * \param[in]   string      the input string
 * \param[out]  b           number of the first non white character
 * \param[out]  e           number of the last non white character
 *----------------------------------------------------------------------------*/

static void
_bounds(const char *str,
       int &start,
       int &end)
{
  int len = strlen(str);
  start = 0;
  end = len - 1;

  while (start <= end && isspace(str[start]))
    start++;
  while (end >= start && isspace(str[end]))
    end--;
}

/*----------------------------------------------------------------------------
 * \brief Check if a line starts with a skip character
 *----------------------------------------------------------------------------*/

static int
_is_skipped(const char* line) {
  int start, end;
  _bounds(line, start, end);
  return (strchr(skip_chars, line[start]) != nullptr);
}

/*----------------------------------------------------------------------------
 * \brief  Find next validated line
 * \param[in]   fp              File
 * \param[in]   current_line    current line
 * \param[out]  source_file     the characters of the line
 *----------------------------------------------------------------------------*/

static int
_find_next_line(FILE        *fp,
                char        *current_line,
                const char  *source_file)
{

  while (fgets(current_line, line_len, fp)) {
    size_t len = strlen(current_line);
    if (len > 0 && current_line[len - 1] == '\n')
      current_line[len - 1] = '\0';

    int start, end;
    _bounds(current_line, start, end);

    if (start <= end && !_is_skipped(current_line))
      return 0; // Valid line found
  }

  if (feof(fp))
    return 1; // end of file
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Unexpected read error on file %s\n"
                "Read error."), source_file);

  return -1;
}

/*----------------------------------------------------------------------------
 * \brief Reads a file having in each significative line a file name
 *  it returns then as 'the_list' the list of lines read
 *  a line is significative if it's first char is not / or # or !
 *  The following 3 lines give an example from which one must remove the
 *  first two characters.
 *  /list of files
 *   profile_one.txt
 *   profile_two.txt
 *
 *  Beware that blank lines or lines starting with blanks
 *  + comment_char areNOT ignored.
 *
 * \param[in]   filename      the name of file with list of file names
 *----------------------------------------------------------------------------*/

static void
_read_files_list(const char *filename)
{
  char line[line_len];

  int count = 0;
  FILE *f = fopen(filename, "r");

  while (_find_next_line(f, line, filename) == 0)
    ++count;

  fclose(f);

  number_of_files = count;

  CS_MALLOC(imbrication_files, number_of_files, char*);

  f = fopen(filename, "r");

  int idx = 0;
  while (_find_next_line(f, line, filename) == 0 && idx < count) {
    char *_line = strdup(line);
    imbrication_files[idx] = _line;
    ++idx;
  }

  fclose(f);
}

/*----------------------------------------------------------------------------
 * Reads a meteo_file for code_saturne Atmospheric Physics option
 *
 * They contain an arbitrary number (>=1) of sections having the following
 * structure.
 * Comment lines start with a slash / as first character
 *
 * yyyy, dd, hh, mm, ss
 * xpos, ypos
 * ground pressure
 * nt (thermal profile dimension)
 * nt lines of
 * zt, temp_c, qw(kg/kg), Ndrops(1/cm3)
 * nd (thermal profile dimension)
 * nd lines of
 * zd, u, v, k, eps
 *
 * WARNINGS:
 * Beware that all dimensions nt,nd must be the same as the first one.
 *
 * Beware that blank lines or lines starting with blanks+ comment_char are
 * NOT ignored.
 *
 * \param[in]   meteo_file      "meteo" file name
 *----------------------------------------------------------------------------*/

static void
_read_meteo_file(const char  *meteo_file)
{
  char current_line[1024];
  int count_lines = 0;
  int first = 0, last = 0;
  int nt = -1, nd = -1, ns = -1;  // ns unused in this extract

  cs_atmo_model_t atmo_model_flag
    = static_cast<cs_atmo_model_t>(cs_glob_physical_model_flag[CS_ATMOSPHERIC]);

  /* Count calls */
  static int file_count = 0;
  file_count++;

  FILE *unilog = fopen(meteo_file, "r");

  /* First loop on file for counting purposes
   * ---------------------------------------- */

  while (true) {
    const int status = _find_next_line(unilog, current_line, meteo_file);
    if (status == 1)
      break;
    count_lines++;
    _bounds(current_line, first, last);
    if (count_lines == 4) {
      int l_iostat = sscanf(current_line, "%d", &nt);
      if (l_iostat != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("Unexpected read error (2) on line: %s\n"
                    "found in file: %s\n"
                    "check and correct your data files\n"
                    "all calculations will be stopped."),
                  current_line, meteo_file);
      if (thermal_profile_dim < 0)
        thermal_profile_dim = nt;
      if (thermal_profile_dim != nt)
        bft_error(__FILE__, __LINE__, 0,
                  _("All files should provide the same nt for thermal profiles\n"
                    "nt read on file >%s< is %d\n"
                    "Global nt read on first file is %d."),
                  meteo_file, nt, thermal_profile_dim);
    }

    if (count_lines > 3) {
      if (count_lines == 4 + nt + 1) {
        if (sscanf(current_line, "%d", &nd) != 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("while reading nd (dynamic profiles dim)\n"
                      "unexpected read error (3) on line: %s\n"
                      "found in file: %s"),
                    current_line, meteo_file);
        if (dynamic_profile_dim < 0)
          dynamic_profile_dim = nd;
        if (dynamic_profile_dim != nd) {
          bft_printf("All files should provide the same nd for dynamic profiles.\n");
          _bounds(meteo_file, first, last);
          bft_error(__FILE__, __LINE__, 0,
                    _("nd read on file >%.*s< is %d\n"
                      "global nd read on first file is %d."),
                    last - first + 1, meteo_file + first,
                    nd, dynamic_profile_dim);
        }
      }
    }

  } // end loop

  if ((count_lines % (nt + nd + 5)) != 0) {
    _bounds(meteo_file, first, last);

    bft_error(__FILE__, __LINE__, 0,
              _("Read_meteo_file encountered an error on file >%.*s<\n"
                "all the sections date, pos, pressure+profiles should"
                "have the same length.\n"
                "Yhis length should be multiple of 5+nt(thermal dim)+nd"
                "(dynamic dim): %d\n"
                "but one found %d lines of data\n"
                "Probable error cause: the actual and expected length of"
                " the profiles differ\n"
                "All thermal profiles should have %d lines.\n"
                "All dynamic profiles should have %d lines."),
              last - first + 1, meteo_file + first,
              5 + nt + nd, count_lines, nt, nd);
  }
  else {
    ns = count_lines / (nt + nd + 5);

    if (sections_per_file < 1)
      sections_per_file = ns;
    if (ns != sections_per_file) {
      _bounds(meteo_file, first, last);
      bft_error(__FILE__, __LINE__, 0,
                _("Read_meteo_file encountered an error on file >%.*s<\n"
                  "all the files should contain the same number of sections"
                  " date, pos, pressure+profiles.\n"
                  "Number of sections of current file is %d.\n"
                  "Number of sections given by first file is %d."),
                last - first + 1, meteo_file + first, ns, sections_per_file);
    }
  }

  /* Here the end of file is reached and nt, nd are consistent with
   * global values as well as the number of sections in the file.
   * --------------------------------------------------------------
   * allocation of arrays for reading data (only the first time)
   * -------------------------------------------------------------- */

  if (file_count == 1)
    _allocate_all();

  fclose(unilog);

  /* -------------------------------------------------------------------
   * second loop on file : one can safely read data to load the profiles
   * ------------------------------------------------------------------- */
  _bounds(meteo_file, first, last);
  if (_atmo_imbrication.imbrication_verbose)
    bft_printf("Reopening the file '%.*s'\n",
               last - first + 1, meteo_file + first);

  unilog = fopen(meteo_file, "r");
  if (!unilog) {
    _bounds(meteo_file, first, last);
    bft_error(__FILE__, __LINE__, 0,
              _("%s: could not open file '%.*s'\n"), __func__,
              last - first + 1, meteo_file + first);
  }

  char line_buf[1024];

  count_lines = 0;
  int section_count = 1;
  int section_length = nt + nd + 5;

  /* loop on file */
  while (true) {
    if (_atmo_imbrication.imbrication_verbose)
      bft_printf("Section count = %d\n"
                 "File count = %d\n", section_count, file_count);
    const int status = _find_next_line(unilog, line_buf, meteo_file);
    if (status == 1) break;
    count_lines++;
    const int mod = count_lines % section_length;
    const int idx = section_count-1+(file_count-1)*sections_per_file;
    if (mod == 1) {
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("Reading years, ord, hours, minutes, seconds in: %s\n",
                   line_buf);
      int sc = sscanf(line_buf, "%d, %d, %d, %d, %lf",
                      &years[idx],
                      &ordinals[idx],
                      &hours[idx],
                      &minutes[idx],
                      &seconds[idx]);
      if (sc != 5)
        bft_error(__FILE__, __LINE__, 0,
                _("While reading years, ord, hours, minutes, seconds\n"
                  "unexpected read error (4) on line: %s\n"
                  "found in file: %s"),
                  line_buf, meteo_file);
    }

    else if (mod == 2) {
     _bounds(line_buf, first, last);
     if (_atmo_imbrication.imbrication_verbose)
       bft_printf("Reading xpos ypos in: %s\n", line_buf);
     int sc = sscanf(line_buf, "%lf %lf",
                     &xpos[idx],
                     &ypos[idx]);
     if (sc != 2)
       bft_error(__FILE__, __LINE__, 0,
                 _("While reading hours, minutes, seconds\n"
                   "unexpected read error (5) on line: %.*s\n"
                   "found in file: %s"),
                 last, line_buf, meteo_file);
    }

    else if (mod == 3) {
      _bounds(line_buf, first, last);
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("reading ground pressure in: %s\n", line_buf);
      const int sc = sscanf(line_buf, "%lf", &ground_pressure[idx]);
      if (sc != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("While reading ground pressure\n"
                    "unexpected read error (6) on line: %.*s\n"
                    "found in file: %s"), last, line_buf, meteo_file);
    }

    for (int i = 1; i <= nt; ++i) {
      if (mod != 4 + i)
        continue;
      const int id = i-1 + (section_count-1)*thermal_profile_dim
                   + (file_count-1)*thermal_profile_dim*sections_per_file;
      _bounds(line_buf, first, last);
      if (atmo_model_flag == CS_ATMO_HUMID) {
        if (_atmo_imbrication.imbrication_verbose)
          bft_printf("Reading zt, tempC, qw Nc in: %.*s\n", last, line_buf);
        int sc = sscanf(line_buf, "%lf %lf %lf %lf",
                        &zt[id],
                        &temp_c[id],
                        &qw[id],
                        &Nc[id]);
        if (sc != 4)
          bft_error(__FILE__, __LINE__, 0,
                    _("While reading ground pressure\n"
                      "unexpected read error (7) on line: %.*s\n"
                      "found in file: %s"), last, line_buf, meteo_file);
      }
      else if (atmo_model_flag < CS_ATMO_HUMID) {
        if (_atmo_imbrication.imbrication_verbose)
          bft_printf("reading zt,tempC,qw in: %s\n", line_buf);
        int sc = sscanf(line_buf, "%lf %lf %lf",
                        &zt[id],
                        &temp_c[id],
                        &qw[id]);

        if (sc != 3)
          bft_error(__FILE__, __LINE__, 0,
                    _("While reading zt, tempC, qw\n"
                      "unexpected read error (8) on line: %.*s\n"
                      "found in file: %s"), last, line_buf, meteo_file);
      }

    }

    for (int i = 1; i <= nd; ++i) {
      const int target = (5 + nt + i) % section_length;
      if (mod != target)
        continue;
      const int id = i-1 + (section_count-1)*dynamic_profile_dim
                   + (file_count-1)*dynamic_profile_dim*sections_per_file;
      _bounds(line_buf, first, last);
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("Reading u,v,tke,eps in: %s\n", line_buf);
      int sc = sscanf(line_buf, "%lf %lf %lf %lf %lf",
                      &zd[id],
                      &vel_u[id],
                      &vel_v[id],
                      &tke[id],
                      &eps[id]);
      if (sc != 5)
        bft_error(__FILE__, __LINE__, 0,
                  _("While reading u, v, tke, eps\n"
                    "unexpected read error (9) on line: %.*s\n"
                    "found in file: %s"), last, line_buf, meteo_file);
    }

    if (mod == 0)
      section_count ++;

  } // end loop

  fclose(unilog);

  if (_atmo_imbrication.imbrication_verbose)
    bft_printf("Read_loop is finished.\n");
}

/*----------------------------------------------------------------------------
 * \brief Checks the time variables to ensure the chronology
 *----------------------------------------------------------------------------*/

static void
_check_chronologies(void)
{
  cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  cs_real_t sim_time;

  /* some checking on chronologies must be done
   * they must be synchronized */
  for (int i = 1; i < number_of_files; i++) {
    for (int j = 0; j < sections_per_file; j++) {
      bool err_chronologies = false;
      const int idx = j + i*sections_per_file;
      if (years[idx] != years[j])
        err_chronologies = true;
      if (ordinals[idx] != ordinals[j])
        err_chronologies = true;
      if (hours[idx] != hours[j])
        err_chronologies = true;
      if (minutes[idx] != minutes[j])
        err_chronologies = true;
      if (fabs(seconds[idx] - seconds[j]) > cs_math_epzero)
        err_chronologies = true;
      if (err_chronologies)
        bft_error(__FILE__, __LINE__, 0,
                  _("the chronologies of the different"
                    " profiles are not synchronized\n"
                    "faulty file: %s\n"
                    "faulty date: %d %d %d %d %lf\n"
                    "should be equal to date: %d %d %d %d %lf\n"
                    "defined in file: %s\n"
                    "section: %d"), imbrication_files[i], years[idx],
                  ordinals[idx], hours[idx], minutes[idx],
                  seconds[idx], years[j], ordinals[j],
                  hours[j], minutes[j], seconds[j],
                  imbrication_files[0], j+1);
    }
  }

  // they must be in the natural temporal order
  const int size = sections_per_file*number_of_files;
  CS_MALLOC(times, size, cs_real_t);

  if (at_opt->syear < 0) {
    at_opt->syear = years[0];
    at_opt->squant  = ordinals[0];
    at_opt->shour = hours[0];
    at_opt->smin  = minutes[0];
    at_opt->ssec  = seconds[0];

    sim_time = _yo2j(years[0], ordinals[0]) * 86400.0
             + hours[0]   * 3600.0
             + minutes[0] * 60.0
             + seconds[0];
  }
  else {
    sim_time = _yo2j(at_opt->syear, at_opt->squant) * 86400.0
             + at_opt->shour * 3600.0
             + at_opt->smin  * 60.0
             + at_opt->ssec;
  }

  for (int i = 0; i < number_of_files; i++) {
    for (int j = 0; j < sections_per_file; j++) {
      const int idx = j + i*sections_per_file;
      times[idx] = _yo2j(years[idx], ordinals[idx]) * 86400.0
                 + hours[idx]   * 3600.0
                 + minutes[idx] * 60.0
                 + seconds[idx];
    }
  }

  for (int i = 0; i < number_of_files; i++)
    for (int j = 0; j < sections_per_file; j++) {
      const int idx = j + i*sections_per_file;
      times[idx] -= sim_time;
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("simulation times: %lf\n", times[idx]);
    }

  for (int i = 0; i < number_of_files; i++) {
    const int ii = i*sections_per_file;
    for (int j = 1; j < sections_per_file; j++) {
      const int idx = j + i*sections_per_file;
      if (times[idx] <= times[ii])
        bft_error(__FILE__, __LINE__, 0,
                  _("the chronologies of the different profiles are not in order\n"
                    "faulty file: %s\n"
                    "faulty date: %d %d %d %d %lf\n"
                    "defined in section %d\n"
                    "should be posterior to date: %d %d %d %d %lf\n"
                    "defined in section %d"), imbrication_files[i], years[idx],
                  ordinals[idx], hours[idx], minutes[idx],
                  seconds[idx], j+1, years[ii], ordinals[ii],
                  hours[ii], minutes[ii], seconds[ii], 1);
    }
  }

  for (int i = 0; i < number_of_files; i++)
    for (int j = 0; j < sections_per_file; j++) {
      const int idx = j + i*sections_per_file;
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("simulation times: %lf\n", times[idx]);
    }
}

/*----------------------------------------------------------------------------
 * \brief Check that the profiles position is the same over time
 *----------------------------------------------------------------------------*/

static void
_check_positions(void)
{
  for (int i = 0; i < number_of_files; i++) {
    for (int j = 1; j < sections_per_file; j++) {
      const int idx = j + i*sections_per_file;
      if (fabs(xpos[idx] - xpos[i*sections_per_file]) > cs_math_epzero)
        bft_error(__FILE__, __LINE__, 0,
                  _("the x-positions of the profiles in file %s\n"
                    "are not consistent (vary in time)\n"
                    "faulty section is : %d\n"
                    " xpos(1)=%lf\n"
                    " xpos(%d)=%lf"), imbrication_files[i],
                  j+1, xpos[i*sections_per_file],
                  j+1, xpos[idx]);

      if (fabs(ypos[idx] - ypos[i*sections_per_file]) > cs_math_epzero)
        bft_error(__FILE__, __LINE__, 0,
                  _("the y-positions of the profiles in file %s\n"
                    "are not consistent (vary in time)\n"
                    "faulty section is : %d\n"
                    " ypos(1)=%lf\n"
                    " ypos(%d)=%lf"), imbrication_files[i],
                  j+1, ypos[i*sections_per_file],  j+1, ypos[idx]);
    }
  }

  for (int i = 0; i < number_of_files; i++)
    for (int k = 0; k < number_of_files; k++)
      if (k != i)
        if (   (     fabs(xpos[i * sections_per_file] - xpos[k * sections_per_file])
                <= cs_math_epzero)
            && (   fabs(ypos[i * sections_per_file] - ypos[k * sections_per_file])
                <= cs_math_epzero))
          bft_error(__FILE__, __LINE__, 0,
                    _("The positions given of some profiles are not consistent\n"
                      "The positions of the profiles in file %s\n"
                      "and the positions of the profiles in file %s\n"
                      "are equal."), imbrication_files[i], imbrication_files[k]);
}

/*----------------------------------------------------------------------------
 * \brief Check that the profiles vertical grids heights
 *        are strictly increasing
 *----------------------------------------------------------------------------*/

static void
_check_altitudes(void)
{
  for (int i = 0; i < number_of_files; i++) {
    for (int j = 0; j < sections_per_file; j++) {
      for (int k = 1; k < thermal_profile_dim; k++) {
        int idx_km1 = (k-1) + j * thermal_profile_dim
          + i * thermal_profile_dim * sections_per_file;
        int idx_k   = k     + j * thermal_profile_dim
          + i * thermal_profile_dim * sections_per_file;

        if (zt[idx_km1] >= zt[idx_k]) {
          for (int kk = 0; kk < thermal_profile_dim; kk++) {
            int idx_kk = kk + j * thermal_profile_dim
              + i * thermal_profile_dim * sections_per_file;
            bft_printf("k=%d zt=%lf\n", kk+1, zt[idx_kk]);
          }
          bft_error(__FILE__, __LINE__, 0,
                    _("the thermal profile in section %d\n"
                      "of the file '%s'\n"
                      "is not strictly increasing\n"
                      "erroneous level %d with zt = %lf"),
                      j+1, imbrication_files[i], k+1, zt[idx_k]);

        }
      }
    }
  }

  for (int i = 0; i < number_of_files; i++) {
    for (int j = 0; j < sections_per_file; j++) {
      for (int k = 1; k < dynamic_profile_dim; k++) {
        int idx_km1 = (k-1) + j * dynamic_profile_dim
          + i * dynamic_profile_dim * sections_per_file;
        int idx_k   =  k    + j * dynamic_profile_dim
          + i * dynamic_profile_dim * sections_per_file;

        if (zd[idx_km1] >= zd[idx_k]) {
          for (int kk = 0; kk < dynamic_profile_dim; kk++) {
            int idx_kk = kk + j * dynamic_profile_dim
              + i * dynamic_profile_dim * sections_per_file;
            bft_printf("k=%d zd=%lf\n", kk+1, zd[idx_kk]);
          }
          bft_error(__FILE__, __LINE__, 0,
                    _("the dynamic profile in section %d\n"
                      "of the file '%s'\n"
                      "is not strictly increasing\n"
                      "erroneous level %d with zd = %lf"),
                    j+1, imbrication_files[i], k+1, zd[idx_k]);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * \brief Compute the hydrostastic pressure by Laplace integration
 *----------------------------------------------------------------------------*/

static void
_hydrostatic_pressure(void)
{
  cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  static int ih2o = 0;
  const cs_real_t rvsra  = cs_glob_fluid_properties->rvsra;
  const cs_real_t gz     = cs_glob_physical_constants->gravity[2];
  const cs_real_t rair   = cs_glob_fluid_properties->r_pg_cnst;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_atmo_model_t atmo_model_flag
    = static_cast<cs_atmo_model_t>(cs_glob_physical_model_flag[CS_ATMOSPHERIC]);

  if (atmo_model_flag == CS_ATMO_HUMID)
    ih2o = 1;

  const int size_pr
    = thermal_profile_dim*sections_per_file*number_of_files;
  if (pressure == nullptr)
    CS_MALLOC(pressure, size_pr, cs_real_t);

  if (at_opt->hydrostatic_pressure_model == 0) {
    for (int i = 0; i < number_of_files; i++) {
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("hydrostatic_pressure::file: %s\n", imbrication_files[i]);

      for (int j = 0; j < sections_per_file; j++) {
        if (_atmo_imbrication.imbrication_verbose)
          bft_printf("hydrostatic_pressure::section: %d\n"
                     "hydrostatic_pressure::thermal_profile_dim=%d\n",
                     j+1, thermal_profile_dim);
        const int id =   j * thermal_profile_dim
                       + i * thermal_profile_dim * sections_per_file;
        pressure[id] = ground_pressure[j + i * sections_per_file];

        for (int k = 1; k < thermal_profile_dim; k++) {
          if (_atmo_imbrication.imbrication_verbose)
            bft_printf("hydrostatic_pressure::k=%d\n", k+1);

          const int idx_k   = k + j * thermal_profile_dim
                              + i * thermal_profile_dim * sections_per_file;
          const int idx_km1 = k-1 + j * thermal_profile_dim
                              + i * thermal_profile_dim * sections_per_file;
          const cs_real_t tmoy = 0.5 * (temp_c[idx_km1] + temp_c[idx_k]) + tkelvi;

          cs_real_t  q0, q1;
          if (atmo_model_flag == CS_ATMO_HUMID) {
            q0 = fmin(qw[idx_km1],
                        cs_air_yw_sat(temp_c[idx_km1], pressure[idx_km1]));
            q1 = fmin(qw[idx_k],
                      cs_air_yw_sat(temp_c[idx_k], pressure[idx_km1]));
          }
          else {
            q0 = qw[idx_km1];
            q1 = qw[idx_k];
          }

          cs_real_t rho_moy
            = rair * (1.0 + (rvsra - 1.0) * (q0 + q1) / 2.0 * ih2o);
          cs_real_t rap = -fabs(gz) * (zt[idx_k] - zt[idx_km1]) / rho_moy / tmoy;

          pressure[idx_k] = pressure[idx_km1] * exp(rap);
        }
      }
    }
  }
  else {
    for (int i = 0; i < number_of_files; i++) {
      for (int j = 0; j < sections_per_file; j++) {
        int idx_top = (thermal_profile_dim-1)
                    + j * thermal_profile_dim
                    + i * thermal_profile_dim*sections_per_file;

        pressure[idx_top] = 101325.0 *
          pow(288.15 / (288.15 - 6.5e-3 * zt[idx_top]),
              -fabs(gz) / rair / 6.5e-3);

        for (int k = thermal_profile_dim-1; k > 0; k--) {
          int idx_k   = k + j * thermal_profile_dim
                      + i * thermal_profile_dim * sections_per_file;
          int idx_km1 = (k-1) + j * thermal_profile_dim
                      + i * thermal_profile_dim * sections_per_file;

          cs_real_t tmoy = 0.5 * (temp_c[idx_km1] + temp_c[idx_k]) + tkelvi;

          cs_real_t q0, q1;
          if (atmo_model_flag == CS_ATMO_HUMID) {
            q0 = fmin(qw[idx_km1],
                      cs_air_yw_sat(temp_c[idx_km1], pressure[idx_k]));
            q1 = fmin(qw[idx_k],
                      cs_air_yw_sat(temp_c[idx_k], pressure[idx_k]));
          }
          else {
            q0 = qw[idx_km1];
            q1 = qw[idx_k];
          }

          cs_real_t rho_moy
            = rair * (1.0 + (rvsra - 1.0) * (q0 + q1) / 2.0 * ih2o);
          cs_real_t rap = fabs(gz) * (zt[idx_k] - zt[idx_km1]) / rho_moy / tmoy;

          pressure[idx_km1] = pressure[idx_k] * exp(rap);
        }
      }
    }
  }

  for (int i = 0; i < number_of_files; i++) {
    if (_atmo_imbrication.imbrication_verbose)
      bft_printf("hydrostatic_pressure::file: %s\n",
                 imbrication_files[i]);

    for (int j = 0; j < sections_per_file; j++) {
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("hydrostatic_pressure::section: %d\n"
                   "hydrostatic_pressure::date: %d %d %d %d %lf\n",
                   j+1, years[j + i * sections_per_file],
                   ordinals[j + i * sections_per_file],
                   hours[j + i * sections_per_file],
                   minutes[j + i * sections_per_file],
                   seconds[j + i * sections_per_file]);

      for (int k = 0; k < thermal_profile_dim; k++)
        if (_atmo_imbrication.imbrication_verbose) {
          int idx = k + j * thermal_profile_dim
            + i * thermal_profile_dim * sections_per_file;
          bft_printf("hydrostatic_pressure::z,t,p: %lf %lf %lf\n",
                     zt[idx], temp_c[idx], pressure[idx]);
        }
    }
  }

}

/*----------------------------------------------------------------------------
 * \brief Compute the potential_temperature_and_density profiles
 *----------------------------------------------------------------------------*/

static void
_potential_temperature_and_density(void)
{
  static int ih2o = 0;
  const cs_real_t ps = cs_glob_atmo_constants->ps;
  const cs_real_t cp0    = cs_glob_fluid_properties->cp0;
  const cs_real_t rvsra = cs_glob_fluid_properties->rvsra;
  const cs_real_t rair = cs_glob_fluid_properties->r_pg_cnst;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
  const cs_real_t cpvcpa = cs_glob_air_props->cp_v / cs_glob_air_props->cp_a;

  cs_atmo_model_t atmo_model_flag
    = static_cast<cs_atmo_model_t>(cs_glob_physical_model_flag[CS_ATMOSPHERIC]);

  if (atmo_model_flag == CS_ATMO_HUMID)
    ih2o = 1;

  const int size
    = thermal_profile_dim*sections_per_file*number_of_files;
  if (theta == nullptr)
    CS_MALLOC(theta, size, cs_real_t);
  if (density == nullptr)
    CS_MALLOC(density, size, cs_real_t);

  for (int i = 0; i < number_of_files; i++) {
    for (int j = 0; j < sections_per_file; j++) {
      for (int k = 0; k < thermal_profile_dim; k++) {

        int idx = k + j * thermal_profile_dim
          + i * thermal_profile_dim * sections_per_file;

        const cs_real_t rhum = rair * (1.0 + (rvsra - 1.0) * qw[idx] * ih2o);

        if (atmo_model_flag == CS_ATMO_CONSTANT_DENSITY) {
          int idx_ref = j * thermal_profile_dim
                      + i * thermal_profile_dim * sections_per_file;
          density[idx] = pressure[idx_ref] / (temp_c[idx] + tkelvi) / rhum;
        }
        else
          density[idx] = pressure[idx] / (temp_c[idx] + tkelvi) / rhum;

        const cs_real_t rscp
          = (rair / cp0) * (1.0 + (rvsra - cpvcpa) * qw[idx] * ih2o);

        theta[idx] = (temp_c[idx] + tkelvi) * pow(ps / pressure[idx], rscp);
      }
    }
  }

  for (int i = 0; i < number_of_files; i++) {
    if (_atmo_imbrication.imbrication_verbose)
      bft_printf("potential_temperature_and_density::file: %s\n",
                 imbrication_files[i]);

    for (int j = 0; j < sections_per_file; j++) {
      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("potential_temperature_and_density::section: %d\n"
                   "potential_temperature_and_density::date: %d %d %d %d %lf\n",
                   j+1, years[j + i * sections_per_file],
                   ordinals[j + i * sections_per_file],
                   hours[j + i * sections_per_file],
                   minutes[j + i * sections_per_file],
                   seconds[j + i * sections_per_file]);

      for (int k = 0; k < thermal_profile_dim; k++) {
        if (_atmo_imbrication.imbrication_verbose) {
          int idx = k + j * thermal_profile_dim
            + i * thermal_profile_dim * sections_per_file;

          bft_printf("z,t,p,potential_temperature,density:::"
                     " %lf %lf %lf %lf %lf\n",
                     zt[idx],
                     temp_c[idx],
                     pressure[idx],
                     theta[idx], density[idx]);
        }
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * \brief Compute radius of influence
 *----------------------------------------------------------------------------*/

static void
_red_tape(void)
{
  const int size_th = 3*thermal_profile_dim*number_of_files;
  const int size_dyn = 3*dynamic_profile_dim*number_of_files;
  const cs_real_t horizontal_influence_radius
    = _atmo_imbrication.horizontal_influence_radius;
  const cs_real_t vertical_influence_radius
    = _atmo_imbrication.vertical_influence_radius;

  CS_MALLOC(coordinates_th, size_th, cs_real_t);
  CS_MALLOC(coordinates_dyn, size_dyn, cs_real_t);
  CS_MALLOC(influence_param_th, size_th, cs_real_t);

  for (int i = 0; i < number_of_files; i++)
    for (int j = 0; j < thermal_profile_dim; j++) {
      const int idx_1 = 0 + j * 3 + i * 3*thermal_profile_dim;
      const int idx_2 = 1 + j * 3 + i * 3*thermal_profile_dim;
      const int idx_3 = 2 + j * 3 + i * 3*thermal_profile_dim;
      influence_param_th[idx_1] = 1.0/horizontal_influence_radius;
      influence_param_th[idx_2] = 1.0/horizontal_influence_radius;
      influence_param_th[idx_3] = 1.0/vertical_influence_radius;
    }

  CS_MALLOC(influence_param_dyn, size_dyn, cs_real_t);

  for (int i = 0; i < number_of_files; i++)
    for (int j = 0; j < dynamic_profile_dim; j++) {
      const int idx_1 = 0 + j * 3 + i * 3 * dynamic_profile_dim;
      const int idx_2 = 1 + j * 3 + i * 3 * dynamic_profile_dim;
      const int idx_3 = 2 + j * 3 + i * 3 * dynamic_profile_dim;
      influence_param_dyn[idx_1] = 1.0/horizontal_influence_radius;
      influence_param_dyn[idx_2] = 1.0/horizontal_influence_radius;
      influence_param_dyn[idx_3] = 1.0/vertical_influence_radius;
    }

}

/*----------------------------------------------------------------------------
 * \brief Fill thermal and dynamic array coordinates
 *----------------------------------------------------------------------------*/

static void
_fill_coordinates(const char       *name,
                  const int        nb_mes,
                  const int        profile_dim,
                  const cs_real_t  z_profile[],
                  const cs_real_t  ti_value[],
                  int              ones[],
                  cs_real_t        coordinates[])

{

  if (_atmo_imbrication.imbrication_verbose)
    bft_printf("nbmes = %d\n", nb_mes);

  if (profile_dim < 1 || number_of_files < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("In %s:"
                "dimensions of time interpolated %s are not consistent\n"
                "expected dimensions are: (1:%d,1:%d)\n"
                "all calculations will be stopped."),
              __func__, name, profile_dim, number_of_files);

  for (int i = 0; i < number_of_files; i++) {
    for (int j = 0; j < profile_dim; j++) {
      const cs_real_t x_i = xpos[i * sections_per_file];
      const cs_real_t y_i = ypos[i * sections_per_file];
      const cs_real_t ti_z = z_profile[j + i * profile_dim];

      const int base = j * 3 + i * 3 * profile_dim;
      coordinates[0 + base] = x_i;
      coordinates[1 + base] = y_i;
      coordinates[2 + base] = ti_z;

      if (_atmo_imbrication.imbrication_verbose)
        bft_printf("j=%d, i=%d, x=%f, y=%f, z=%f, %s=%f\n",
                   j+1, i+1,
                   coordinates[0 + base],
                   coordinates[1 + base],
                   coordinates[2 + base],
                   name,
                   ti_value[j + i*profile_dim]);

      ones[j + i * profile_dim] = 1;

    }
  }
}

/*----------------------------------------------------------------------------
 * \brief Fill in a measures set structure with an array of measures
 *----------------------------------------------------------------------------*/

static void
_measures_set_map_values(const char       *name,
                         const int        nb_mes,
                         const int        profile_dim,
                         const int        imbrication_id,
                         const cs_real_t  z_profile[],
                         const cs_real_t  ti_value[],
                         const cs_real_t  influence_param[],
                         cs_real_t        coordinates[])
{
  int *ones = nullptr;
  CS_MALLOC(ones, nb_mes, int);

  _fill_coordinates(name,
                    nb_mes,
                    profile_dim,
                    z_profile,
                    ti_value,
                    ones,
                    coordinates);

  cs_measures_set_t *ms = cs_measures_set_by_id(imbrication_id);
  cs_measures_set_map_values(ms,
                             nb_mes,
                             ones,
                             ones,
                             coordinates,
                             ti_value,
                             influence_param);

  CS_FREE(ones);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * \brief Prepare data for imbrication by reading meteo files
 *
 * Warning : the list of files is supposed to be "imbrication_files_list.txt"
 *----------------------------------------------------------------------------*/

void
cs_activate_imbrication(void)
{
  cs_atmo_model_t atmo_model_flag
    = static_cast<cs_atmo_model_t>(cs_glob_physical_model_flag[CS_ATMOSPHERIC]);

  cs_log_separator(CS_LOG_DEFAULT);
  cs_log_printf(CS_LOG_DEFAULT,
                _("Atmospheric Imbrication:\n"));
  cs_log_separator(CS_LOG_DEFAULT);

  imbrication_files_list = "imbrication_files_list.txt";
  _read_files_list(imbrication_files_list);
  bft_printf("number_of_files            : %d\n", number_of_files);
  for (int ii = 0; ii < number_of_files; ii++) {
    const char *_imbrication_files = (char *)imbrication_files[ii];
    int first, last;
    _bounds(_imbrication_files, first, last);
    if (_atmo_imbrication.imbrication_verbose)
      bft_printf("file number %s\n", imbrication_files[ii]);
    _read_meteo_file(imbrication_files[ii]);
  }

  if (_atmo_imbrication.imbrication_verbose) {
    for (int ii = 0; ii < number_of_files; ii++) {
      const char *_imbrication_files = (char *)imbrication_files[ii];
      int first, last;
      _bounds(_imbrication_files, first, last);
      bft_printf("file number %d = '%.*s'\n",
                 ii+1, last-first+1, _imbrication_files+first);
      bft_printf("number of sections per file: %d\n", sections_per_file);
      for (int jj = 0; jj < sections_per_file; jj++) {
        const int idx = jj + ii*sections_per_file;
        bft_printf("date : %d, %d, %d, %d, %lf\n", years[idx],ordinals[idx],
                   hours[idx],minutes[idx],seconds[idx]);
        bft_printf("xpos, ypos : %lf, %lf\n", xpos[idx], ypos[idx]);
        bft_printf("ground_pressure : %lf\n", ground_pressure[idx]);
        bft_printf("thermal profiles dim : %d\n", thermal_profile_dim);
        for (int kk = 0; kk < thermal_profile_dim; kk++) {
          const int id
            = kk + jj*thermal_profile_dim + ii*thermal_profile_dim*sections_per_file;
          if (atmo_model_flag == CS_ATMO_HUMID)
            bft_printf("z, temp, qw, nc = %lf %lf %lf %lf\n", zt[id],
                       temp_c[id],qw[id],Nc[id]);
          else if (atmo_model_flag == CS_ATMO_DRY)
            bft_printf("z, temp, qw = %lf %lf %lf\n", zt[id],
                       temp_c[id],qw[id]);
          else
            bft_printf("z, temp = %lf %lf\n", zt[id], temp_c[id]);
        } // loop_on_thermal_profile
        bft_printf("     dynamic profiles dim: %d\n", dynamic_profile_dim);
        for (int kk = 0; kk < dynamic_profile_dim; kk++) {
          const int id
            = kk + jj*dynamic_profile_dim + ii*dynamic_profile_dim*sections_per_file;
         bft_printf("z, u, v, k, eps = %lf %lf %lf %lf %lf\n", zd[id], vel_u[id],
                    vel_v[id], tke[id], eps[id]);
        } // loop_on_dynamic_profile
      } // loop_on_sections
    } // loop on file
  } // imbrication_verbose

  /* Some checking on chronologies must be done
   * ------------------------------------------ */
  _check_chronologies();

  /* Some checking on positions of the profiles
   * ------------------------------------------ */
  _check_positions();
  _check_altitudes();

  /* Reading terminated: some calculations are done
   * ------------------------------------------------
   * calculating pressure by hydrostatic (Laplace) integration
   * ------------------------------------------------ */
  _hydrostatic_pressure();

   /* Calculating potential temperature and density
    * --------------------------------------------- */
  _potential_temperature_and_density();
}

/*----------------------------------------------------------------------------
 * \brief Prepare for the cressman interpolation of the variables
 *
 * \param[in]  the_time        current time
 *----------------------------------------------------------------------------*/

void
cs_summon_cressman(cs_real_t the_time)
{
  cs_atmo_model_t atmo_model_flag
    = static_cast<cs_atmo_model_t>(cs_glob_physical_model_flag[CS_ATMOSPHERIC]);

  static bool first_call = true;

  cs_measures_set_t *ms = nullptr;

  if (first_call) {
    if (_atmo_imbrication.cressman_u) {
      ms = cs_measures_set_create("u", 0, 1, false);
      _atmo_imbrication.id_u = ms->id;
    }
    if (_atmo_imbrication.cressman_v) {
      ms = cs_measures_set_create("v", 0, 1, false);
      _atmo_imbrication.id_v = ms->id;
    }
    if (_atmo_imbrication.cressman_tke) {
      ms = cs_measures_set_create("tke", 0, 1, false);
      _atmo_imbrication.id_tke = ms->id;
    }
    if (_atmo_imbrication.cressman_eps) {
      ms = cs_measures_set_create("eps", 0, 1, false);
      _atmo_imbrication.id_eps = ms->id;
    }
    if (_atmo_imbrication.cressman_theta
        && atmo_model_flag >= CS_ATMO_DRY) {
      ms = cs_measures_set_create("theta", 0, 1, false);
      _atmo_imbrication.id_theta = ms->id;
    }
    if (   _atmo_imbrication.cressman_qw
        && atmo_model_flag >= CS_ATMO_HUMID) {
      ms = cs_measures_set_create("qw", 0, 1, false);
      _atmo_imbrication.id_qw = ms->id;
    }
    if (   _atmo_imbrication.cressman_nc
        && atmo_model_flag >= CS_ATMO_HUMID) {
      ms = cs_measures_set_create("nc", 0, 1, false);
      _atmo_imbrication.id_nc = ms->id;
    }

    _red_tape();
    first_call = false;

  } // first call

  _interpolate_all_profiles(the_time);

  if (_atmo_imbrication.cressman_u) {
    const int nbmes = dynamic_profile_dim * number_of_files;

    _measures_set_map_values("u",
                             nbmes,
                             dynamic_profile_dim,
                             _atmo_imbrication.id_u,
                             ti_zd,
                             ti_vel_u,
                             influence_param_dyn,
                             coordinates_dyn);
  }

  if (_atmo_imbrication.cressman_v) {
    const int nbmes = dynamic_profile_dim * number_of_files;
    _measures_set_map_values("v",
                             nbmes,
                             dynamic_profile_dim,
                             _atmo_imbrication.id_v,
                             ti_zd,
                             ti_vel_v,
                             influence_param_dyn,
                             coordinates_dyn);
  }

  if (_atmo_imbrication.cressman_tke) {
    const int nbmes = dynamic_profile_dim * number_of_files;
    _measures_set_map_values("tke",
                             nbmes,
                             dynamic_profile_dim,
                             _atmo_imbrication.id_tke,
                             ti_zd,
                             ti_tke,
                             influence_param_dyn,
                             coordinates_dyn);

  }

  if (_atmo_imbrication.cressman_eps) {
    const int nbmes = dynamic_profile_dim * number_of_files;
    _measures_set_map_values("eps",
                             nbmes,
                             dynamic_profile_dim,
                             _atmo_imbrication.id_eps,
                             ti_zd,
                             ti_eps,
                             influence_param_dyn,
                             coordinates_dyn);
  }

  if (   _atmo_imbrication.cressman_theta
      && (   atmo_model_flag == CS_ATMO_DRY
          || atmo_model_flag == CS_ATMO_HUMID)) {
    const int nbmes = thermal_profile_dim * number_of_files;
    _measures_set_map_values("theta",
                             nbmes,
                             thermal_profile_dim,
                             _atmo_imbrication.id_theta,
                             ti_zt,
                             ti_theta,
                             influence_param_th,
                             coordinates_th);
  }

  if (   _atmo_imbrication.cressman_qw
      && atmo_model_flag == CS_ATMO_HUMID) {
    const int nbmes = thermal_profile_dim * number_of_files;
    _measures_set_map_values("qw",
                             nbmes,
                             thermal_profile_dim,
                             _atmo_imbrication.id_qw,
                             ti_zt,
                             ti_qw,
                             influence_param_th,
                             coordinates_th);

  }

  if (   _atmo_imbrication.cressman_nc
      && atmo_model_flag == CS_ATMO_HUMID) {
    const int nbmes = thermal_profile_dim * number_of_files;
    _measures_set_map_values("nc",
                             nbmes,
                             thermal_profile_dim,
                             _atmo_imbrication.id_nc,
                             ti_zt,
                             ti_Nc,
                             influence_param_th,
                             coordinates_th);
  }
}

/*----------------------------------------------------------------------------
 * \brief  Final step for free arrays imbrication
 *----------------------------------------------------------------------------*/

void
cs_finalize_imbrication(void)
{
  cs_measures_sets_destroy();

  CS_FREE(xpos);
  CS_FREE(ypos);
  CS_FREE(years);
  CS_FREE(hours);
  CS_FREE(minutes);
  CS_FREE(seconds);
  CS_FREE(ti_theta);
  CS_FREE(ti_density);
  CS_FREE(ground_pressure);
  CS_FREE(zt);
  CS_FREE(qw);
  CS_FREE(Nc);
  CS_FREE(zd);
  CS_FREE(tke);
  CS_FREE(eps);
  CS_FREE(vel_u);
  CS_FREE(vel_v);
  CS_FREE(temp_c);
  CS_FREE(times);
  CS_FREE(pressure);
  CS_FREE(theta);
  CS_FREE(density);
  CS_FREE(ti_zt);
  CS_FREE(ti_temp_c);
  CS_FREE(ti_qw);
  CS_FREE(ti_Nc);
  CS_FREE(ti_zd);
  CS_FREE(ti_vel_u);
  CS_FREE(ti_vel_v);
  CS_FREE(ti_tke);
  CS_FREE(ti_eps);
  CS_FREE(ordinals);
  CS_FREE(ti_pressure);
  CS_FREE(coordinates_th);
  CS_FREE(coordinates_dyn);
  CS_FREE(imbrication_files);
  CS_FREE(influence_param_th);
  CS_FREE(influence_param_dyn);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
