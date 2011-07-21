/*============================================================================
 * Definitions, global variables, and base functions
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2008  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

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
