/*============================================================================
 * Definitions, global variables, and base functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_config_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Names of (multiple) "nodal" element types */

const char  *fvm_elements_type_name[] = {N_("edges"),
                                         N_("triangles"),
                                         N_("quadrangles"),
                                         N_("simple polygons"),
                                         N_("tetrahedra"),
                                         N_("pyramids"),
                                         N_("prisms"),
                                         N_("hexahedra"),
                                         N_("simple polyhedra")};

/* Names of (single) "nodal" element types */

const char  *fvm_element_type_name[] = {N_("edge"),
                                        N_("triangle"),
                                        N_("quadrangle"),
                                        N_("simple polygon"),
                                        N_("tetrahedron"),
                                        N_("pyramid"),
                                        N_("prism"),
                                        N_("hexahedron"),
                                        N_("simple polyhedron")};

/* Sizes associated with datatypes */

const size_t  fvm_datatype_size[] = {0,
                                     1,
                                     sizeof(float),
                                     sizeof(double),
                                     4,
                                     8,
                                     4,
                                     8};

const char   *fvm_datatype_name[] = {"",
                                     "char",
                                     "float",
                                     "double",
                                     "int32",
                                     "int64",
                                     "uint32",
                                     "uint64"};

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
