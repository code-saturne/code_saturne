/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
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

#ifndef __CS_OPTS_H__
#define __CS_OPTS_H__

/*============================================================================
 * Parsing of program arguments and associated initializations
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Code_Saturne base options structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char          *app_name;      /* Application name */

  /* Redirection of standard output */

  cs_int_t       ilisr0;        /* Redirection for rank 0
                                   (0: not redirected;
                                   1: redirected to "listing" file) */
  cs_int_t       ilisrp;        /* Redirection for ranks > 0
                                   (0: not redirected;
                                   1: redirected to "listing_n*" file;
                                   2: redirected to "/dev/null", suppressed) */

  /* MPI-IO mode */

  int            mpi_io_mode;   /* MPI-IO mode:
                                   0: no MPI-IO
                                   1: MPI-IO with explicit offsets
                                   2: MPI-IO with individual file pointers */

  /* Other options */

  cs_bool_t      preprocess;    /* Mesh preprocessing mode */
  cs_bool_t      verif;         /* Mesh quality verification mode */

  int            benchmark;     /* Benchmark mode:
                                   0: not used;
                                   1: timing (CPU + Walltime) mode
                                   2: MPI trace-friendly mode */

  /* Server socket for SYRTHES 3 coupling */

  int           syr_socket;     /* socket number to use, or -1 */

  /* Connection with YACS or proxy */

  char          *yacs_module;   /* Path to YACS module */

  char          *proxy_socket;  /* Name of proxy socket */
  int            proxy_key;     /* Key for connection to proxy */

} cs_opts_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define options and call some associated initializations
 * based on command line arguments
 *
 * parameters:
 *   argc  --> number of command line arguments
 *   argv  --> array of command line arguments
 *   opts  <-- options structure
 *----------------------------------------------------------------------------*/

void
cs_opts_define(int         argc,
               char       *argv[],
               cs_opts_t  *opts);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_OPTS_H__ */
