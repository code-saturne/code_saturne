/*============================================================================
 * Definitions, global variables, and base functions
 *============================================================================*/

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

#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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

/*----------------------------------------------------------------------------*/

END_C_DECLS
