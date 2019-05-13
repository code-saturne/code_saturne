#ifndef __FVM_TO_VTK_HISTOGRAM_H__
#define __FVM_TO_VTK_HISTOGRAM_H__

/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to png histogram files
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_writer.h"
#include "fvm_to_histogram.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void
fvm_to_vtk_display_histogram_png(cs_real_t                   var_min,
                                 cs_real_t                   var_max,
                                 cs_gnum_t                   count[],
                                 fvm_to_histogram_writer_t  *w,
                                 char                       *var_name);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_TO_VTK_HISTOGRAM_H__ */
