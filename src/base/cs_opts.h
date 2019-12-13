#ifndef __CS_OPTS_H__
#define __CS_OPTS_H__

/*============================================================================
 * Parsing of program arguments and associated initializations
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

  bool           trace;         /* Trace progress to stdout */
  bool           logrp;         /* Redirection for ranks > 0
                                   (false: to "/dev/null", suppressed)
                                    true: to "run_solver_n*.log" file) */

  /* Signal handling */

  bool           sig_defaults;  /* use default signal handlers */

  /* Other options */

  bool           preprocess;    /* Mesh preprocessing mode */
  bool           verif;         /* Mesh quality verification mode */
  int            benchmark;     /* Benchmark mode:
                                   0: not used;
                                   1: timing (CPU + Walltime) mode
                                   2: MPI trace-friendly mode */

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
