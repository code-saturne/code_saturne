/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_DLOPEN)
#include <dlfcn.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_at_plugin.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_at_plugin.c

  \brief Plugin to dynamically load(*.so) the aerosol model (SIREAM)
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Function Definitions
 *============================================================================*/

#if defined(HAVE_DLOPEN)

/*----------------------------------------------------------------------------
 * Get a shared library function pointer
 *
 * parameters:
 *   handle           <-- pointer to shared library (result of dlopen)
 *   name             <-- name of function symbol in library
 *   errors_are_fatal <-- abort if true, silently ignore if false
 *
 * returns:
 *   pointer to function in shared library
 *----------------------------------------------------------------------------*/

static void *
_get_dl_function_pointer(void        *handle,
                         const char  *lib_path,
                         const char  *name,
                         bool         errors_are_fatal)
{
  void  *retval = NULL;
  char  *error = NULL;
  char  *name_ = NULL;

  dlerror();    /* Clear any existing error */

  retval = dlsym(handle, name);
  error = dlerror();

  if (error != NULL) { /* Try different symbol names */
    dlerror();    /* Clear any existing error */
    int _size_ = strlen(name) + strlen("_");
    BFT_MALLOC(name_, _size_ + 1, char);
    strcpy(name_, name);
    strcat(name_, "_");
    retval = dlsym(handle, name_);
    error = dlerror();
    BFT_FREE(name_);
  }

  if (error != NULL && errors_are_fatal)
    bft_error(__FILE__, __LINE__, 0,
              _("Error while trying to find symbol %s in lib %s: %s\n"),
              name,
              lib_path,
              dlerror());

  return retval;
}

#endif /* defined(HAVE_DLOPEN)*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/* Plug-in to get aerosol function
   from SIREAM library (ENPC - INRIA - EDF R&D) */

void CS_PROCF(plug_aerosol, PLUG_AEROSOL)
(
 cs_int_t   *nx,
 cs_int_t   *ny,
 cs_int_t   *nz,
 cs_int_t   *ns,
 cs_real_t  *ts,
 cs_real_t  *dlhumid,
 cs_real_t  *dltemp,
 cs_real_t  *dlpress,
 cs_real_t  *delta_t,
 cs_real_t  *dlconc,
 cs_int_t   *noptions_aer,
 cs_int_t   *option_aer,
 cs_int_t   *ns_aer,
 cs_int_t   *nbin_aer,
 cs_int_t   *ncycle_aer,
 cs_real_t  *bin_bound_aer,
 cs_real_t  *fixed_density_aer,
 cs_real_t  *density_aer,
 cs_int_t   *couples_coag,
 cs_int_t   *first_index_coag,
 cs_int_t   *second_index_coag,
 cs_real_t  *coefficient_coag,
 cs_real_t  *dlconc_aer,
 cs_real_t  *dlnum_aer
)
{

#if defined(HAVE_DLOPEN)
  typedef void (*aerosol_t)(cs_int_t*, cs_int_t*, cs_int_t*, cs_int_t*,
                            cs_real_t*, cs_real_t*, cs_real_t*, cs_real_t*,
                            cs_real_t*, cs_real_t*, cs_int_t*, cs_int_t*,
                            cs_int_t*, cs_int_t*, cs_int_t*, cs_real_t*,
                            cs_real_t*, cs_real_t*, cs_int_t*, cs_int_t*,
                            cs_int_t*, cs_real_t*, cs_real_t*, cs_real_t*);

  void *handle;
  const char lib_path[] = "libsiream.so";

  handle = dlopen(lib_path, RTLD_LAZY);

  bft_error(__FILE__, __LINE__, 0,
            _("Error loading %s: %s."), lib_path, dlerror());

  aerosol_t aerosol = (aerosol_t) _get_dl_function_pointer(handle,
                                                           lib_path,
                                                           "aerosol",
                                                           true);

  aerosol(nx, ny, nz, ns, ts, dlhumid, dltemp, dlpress, delta_t,
          dlconc, noptions_aer, option_aer, ns_aer, nbin_aer, ncycle_aer,
          bin_bound_aer, fixed_density_aer, density_aer, couples_coag,
          first_index_coag, second_index_coag, coefficient_coag, dlconc_aer,
          dlnum_aer);

  dlclose(handle);

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to load: %s\n"), "libsiream.so");

#endif

}

/* Plug-in to get compute_coagulation_coefficient function
   from SIREAM library (ENPC - INRIA - EDF R&D) */

void CS_PROCF(plug_compute_coagulation_coefficient,
              PLUG_COMPUTE_COAGULATION_COEFFICIENT)
(
 cs_int_t   *nbin_aer,
 cs_real_t  *bin_bound,
 cs_int_t   *couple,
 cs_int_t   *first_index,
 cs_int_t   *second_index,
 cs_real_t  *partition_coefficient
)
{
#if defined(HAVE_DLOPEN)
  typedef void (*compute_coagulation_coefficient_t)(cs_int_t*, cs_real_t*,
                                                    cs_int_t*, cs_int_t*,
                                                    cs_int_t*, cs_real_t*);
  void *handle;
  const char lib_path[] = "libsiream.so";

  handle = dlopen(lib_path, RTLD_LAZY);

  bft_error(__FILE__, __LINE__, 0,
            _("Error loading %s: %s."), lib_path, dlerror());

  compute_coagulation_coefficient_t compute_coagulation_coefficient =
  (compute_coagulation_coefficient_t) _get_dl_function_pointer(handle,
                                                               lib_path,
                                      "compute_coagulation_coefficient",
                                                               true);

  compute_coagulation_coefficient(nbin_aer, bin_bound, couple, first_index,
                                  second_index, partition_coefficient);

  dlclose(handle);

#else

  bft_error(__FILE__, __LINE__, 0,
            _("Shared library support not available.\n"
              "Unable to load: %s\n"), "libsiream.so");

#endif

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
