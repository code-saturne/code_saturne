//============================================================================
//
//     This file is part of the Code_Saturne CFD tool.
//
//     Copyright (C) 2006-2011 EDF S.A., France
//
//     contact: saturne-support@edf.fr
//
//     The Code_Saturne CFD tool is free software; you can redistribute it
//     and/or modify it under the terms of the GNU General Public License
//     as published by the Free Software Foundation; either version 2 of
//     the License, or (at your option) any later version.
//
//     The Code_Saturne CFD tool is distributed in the hope that it will be
//     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with the Code_Saturne Kernel; if not, write to the
//     Free Software Foundation, Inc.,
//     51 Franklin St, Fifth Floor,
//     Boston, MA  02110-1301  USA
//
//============================================================================

//============================================================================
// Unit tests.
//============================================================================

#include "cs_config.h"

// System headers

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <calcium.h>

// Library headers

#include "cfd_proxy_defs.h"
#include "cfd_proxy_child.h"
#include "cfd_proxy_comm.h"
#include "cfd_proxy_api.h"

#include "cfd_proxy_loader.h"

//============================================================================
//                      Local structure definitions
//============================================================================

double (* _ddot)(const int N, const double *X, const int incX,
                 const double *Y, const int incY);

int (* cs_main)(int argc, char **argv);

//============================================================================
//                       Dummy Calcium functions
//============================================================================

// Connection functions

int
cp_cd(void  *component,
      char  *s)
{
  printf("cp_cd(%p, s) called\n", component);

  strcpy(s, "Dummy calcium unit test component");
  return 0;
}

int
cp_fin(void  *component,
       int    cont)
{
  printf("cp_fin(%p, %d) called\n", component, cont);

  return 0;
}

// Blocking read functions

int
cp_len(void  *component,
       int    time_dep,
       float *min_time,
       float *max_time,
       int   *iteration,
       char  *var_name,
       int    n_max_vals,
       int   *n_vals,
       int   *vals)
{
  printf("cp_len(%p, %d, %f, %f, %d, %s, %d, n_vals, vals) called\n",
         component, time_dep, *min_time, *max_time, *iteration,
         var_name, n_max_vals);

  *min_time += 0.5;
  *iteration += 1;
  *n_vals = n_max_vals - 1;

  for (int i = 0; i < n_max_vals; i++)
    vals[i] = i+1;

  return 0;
}

int
cp_lre(void   *component,
       int     time_dep,
       float  *min_time,
       float  *max_time,
       int    *iteration,
       char   *var_name,
       int     n_max_vals,
       int    *n_vals,
       float  *vals)
{
  printf("cp_lre(%p, %d, %f, %f, %d, %s, %d, n_vals, vals) called\n",
         component, time_dep, *min_time, *max_time, *iteration,
         var_name, n_max_vals);

  *min_time += 0.5;
  *iteration += 1;
  *n_vals = n_max_vals - 1;

  for (int i = 0; i < n_max_vals; i++)
    vals[i] = i+1;

  return 0;
}

int
cp_ldb(void    *component,
       int      time_dep,
       double  *min_time,
       double  *max_time,
       int     *iteration,
       char    *var_name,
       int      n_max_vals,
       int     *n_vals,
       double  *vals)
{
  printf("cp_ldb(%p, %d, %f, %f, %d, %s, %d, n_vals, vals) called\n",
         component, time_dep, *min_time, *max_time, *iteration,
         var_name, n_max_vals);

  *min_time += 0.5;
  *iteration += 1;
  *n_vals = n_max_vals - 1;

  for (int i = 0; i < n_max_vals; i++)
    vals[i] = i+1;

  return 0;
}

// Write functions

int
cp_een(void  *component,
       int    time_dep,
       float  cur_time,
       int    iteration,
       char  *var_name,
       int    n_vals,
       int   *vals)
{
  printf("cp_een(%p, %d, %f, %d, %s, %d, vals) called\n",
         component, time_dep, cur_time, iteration, var_name, n_vals);

  for (int i = 0; i < n_vals; i++)
    printf("  vals[%d] = %d\n", i+1, vals[i]);

  return 0;
}

int
cp_ere(void   *component,
       int     time_dep,
       float   cur_time,
       int     iteration,
       char   *var_name,
       int     n_vals,
       float  *vals)
{
  printf("cp_ere(%p, %d, %f, %d, %s, %d, vals) called\n",
         component, time_dep, cur_time, iteration, var_name, n_vals);

  for (int i = 0; i < n_vals; i++)
    printf("  vals[%d] = %f\n", i+1, (double)(vals[i]));

  return 0;
}

int
cp_edb(void    *component,
       int      time_dep,
       double   cur_time,
       int      iteration,
       char    *var_name,
       int      n_vals,
       double  *vals)
{
  printf("cp_edb(%p, %d, %f, %d, %s, %d, vals) called\n",
         component, time_dep, cur_time, iteration, var_name, n_vals);

  for (int i = 0; i < n_vals; i++)
    printf("  vals[%d] = %f\n", i+1, vals[i]);

  return 0;
}

//============================================================================
//                      Private Function Definitions
//============================================================================

//============================================================================
//                      Public Function Definitions
//============================================================================

//----------------------------------------------------------------------------
// Main program
//----------------------------------------------------------------------------

int main
(
 int argc,          // Nombre d'arguments dans la ligne de commandes
 char *argv[]       // Tableau des arguments de la ligne de commandes
)
{
  int i;
  int retval = 0;

  // Default definitions of coupling options
  //----------------------------------------

  // Spawn a child process and establish a connection

  for (i = 0; i < argc; i++) {

    if (strcmp(argv[i], "--dir") == 0)
      cfd_proxy_set_dir(argv[i+1]);

    if (strcmp(argv[i], "--exe") == 0) {

      const char *child_argv[] = {"-param", "param.xml"};

      // char *launcher_args[] = {"/usr/bin/mpiexec", "-n", "1"};
      // cfd_proxy_set_launcher(3, launcher_args);

      cfd_proxy_set_exe(argv[i+1]);

      cfd_proxy_set_args_by_list(2, child_argv);

      cfd_proxy_set_args("");

      retval = cfd_proxy_run_all();

      cfd_proxy_printf("returned from spawn (%d)\n", retval);

    }

  }

#if defined(HAVE_DLOPEN)

  for (i = 0; i < argc; i++) {

    if (strcmp(argv[i], "--lib") == 0) {

      const char *child_argv[] = {"-param", "param.xml"};

      cfd_proxy_set_lib(argv[i+1]);

      cfd_proxy_set_args_by_list(2, child_argv);

      cfd_proxy_set_args("-q");

      retval = cfd_proxy_run_all();

    }

  }

#endif

  return retval;
}
