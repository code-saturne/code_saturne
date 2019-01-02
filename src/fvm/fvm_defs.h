#ifndef __FVM_DEFS_H__
#define __FVM_DEFS_H__

/*============================================================================
 * Definitions, global variables, and base functions
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

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Element types
 *----------------------------------------------------------------------------*/

typedef enum {

  FVM_EDGE,               /* Edge */
  FVM_FACE_TRIA,          /* Triangle */
  FVM_FACE_QUAD,          /* Quadrangle */
  FVM_FACE_POLY,          /* Simple Polygon */
  FVM_CELL_TETRA,         /* Tetrahedron */
  FVM_CELL_PYRAM,         /* Pyramid */
  FVM_CELL_PRISM,         /* Prism (pentahedron) */
  FVM_CELL_HEXA,          /* Hexahedron (brick) */
  FVM_CELL_POLY,          /* Simple Polyhedron (convex or quasi-convex) */
  FVM_N_ELEMENT_TYPES     /* Number of element types */

} fvm_element_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Names of (multiple) element types */

extern const char  *fvm_elements_type_name[];

/* Names of (single) element types */

extern const char  *fvm_element_type_name[];

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_DEFS_H__ */
