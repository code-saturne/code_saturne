/*============================================================================
 * Append sections to a nodal representation associated with a mesh
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_nodal_append.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create new section, transferring ownership of the given connectivity
 * and optional parent number arrays to that section.
 *
 * parameters:
 *   n_elements         <-- number of elements in section
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

static fvm_nodal_section_t *
_transfer_to_section(cs_lnum_t       n_elements,
                     fvm_element_t   type,
                     cs_lnum_t       face_index[],
                     cs_lnum_t       face_num[],
                     cs_lnum_t       vertex_index[],
                     cs_lnum_t       vertex_num[],
                     cs_lnum_t       parent_element_num[])
{
  fvm_nodal_section_t  *this_section = NULL;

  this_section = fvm_nodal_section_create(type);

  this_section->n_elements = n_elements;

  /* Connectivity */

  if (type == FVM_CELL_POLY) {
    this_section->_face_index = face_index;
    this_section->_face_num = face_num;
  }

  if (type == FVM_FACE_POLY || type == FVM_CELL_POLY)
    this_section->_vertex_index = vertex_index;

  this_section->_vertex_num = vertex_num;

  this_section->_parent_element_num = parent_element_num;

  /* Shared arrays */

  this_section->face_index = this_section->_face_index;
  this_section->face_num = this_section->_face_num;
  this_section->vertex_index = this_section->_vertex_index;
  this_section->vertex_num = this_section->_vertex_num;
  this_section->parent_element_num = this_section->_parent_element_num;

  /* Connectivity size */

  if (this_section->stride != 0)
    this_section->connectivity_size
      = this_section->n_elements * this_section->stride;

  else if (this_section->type == FVM_FACE_POLY)
    this_section->connectivity_size
      = this_section->vertex_index[this_section->n_elements];

  else if (this_section->type == FVM_CELL_POLY) {
    cs_lnum_t i, _face_num;
    for (i = 0;
         i < this_section->face_index[this_section->n_elements];
         i++) {
      _face_num = CS_ABS(this_section->face_num[i]);
      if (_face_num > this_section->n_faces)
        this_section->n_faces = _face_num;
    }
    this_section->connectivity_size
      = this_section->vertex_index[this_section->n_faces];
  }

  return this_section;
}

/*----------------------------------------------------------------------------
 * Create new section, mapping the given connectivity and optional
 * parent number arrays to that section.
 *
 * parameters:
 *   n_elements         <-- number of elements in section
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

static fvm_nodal_section_t *
_map_to_section(cs_lnum_t       n_elements,
                fvm_element_t   type,
                cs_lnum_t       face_index[],
                cs_lnum_t       face_num[],
                cs_lnum_t       vertex_index[],
                cs_lnum_t       vertex_num[],
                cs_lnum_t       parent_element_num[])
{
  fvm_nodal_section_t  *this_section = NULL;

  this_section = fvm_nodal_section_create(type);

  this_section->n_elements = n_elements;

  /* Connectivity */

  if (type == FVM_CELL_POLY) {
    this_section->face_index = face_index;
    this_section->face_num = face_num;
  }

  if (type == FVM_FACE_POLY || type == FVM_CELL_POLY)
    this_section->vertex_index = vertex_index;

  this_section->vertex_num = vertex_num;

  this_section->parent_element_num = parent_element_num;

  /* Connectivity size */

  if (this_section->stride != 0)
    this_section->connectivity_size
      = this_section->n_elements * this_section->stride;

  else if (this_section->type == FVM_FACE_POLY)
    this_section->connectivity_size
      = this_section->vertex_index[this_section->n_elements];

  else if (this_section->type == FVM_CELL_POLY) {
    cs_lnum_t i, _face_num;
    for (i = 0;
         i < this_section->face_index[this_section->n_elements];
         i++) {
      _face_num = CS_ABS(this_section->face_num[i]);
      if (_face_num > this_section->n_faces)
        this_section->n_faces = _face_num;
    }
    this_section->connectivity_size
      = this_section->vertex_index[this_section->n_faces];
  }

  return this_section;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Append a new section to an existing fvm_nodal mesh, and transfer
 * ownership of the given connectivity and optional parent number arrays to
 * that section.
 *
 * parameters:
 *   this_nodal         <-> nodal mesh structure
 *   n_elements         <-- number of elements to add
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

void
fvm_nodal_append_by_transfer(fvm_nodal_t    *this_nodal,
                             cs_lnum_t       n_elements,
                             fvm_element_t   type,
                             cs_lnum_t       face_index[],
                             cs_lnum_t       face_num[],
                             cs_lnum_t       vertex_index[],
                             cs_lnum_t       vertex_num[],
                             cs_lnum_t       parent_element_num[])
{
  fvm_nodal_section_t  *new_section = NULL;
  int  n_sections = 0;

  assert(this_nodal != NULL);

  n_sections = this_nodal->n_sections;

  /* Create new section */

  BFT_REALLOC(this_nodal->sections, n_sections + 1, fvm_nodal_section_t *);

  new_section = _transfer_to_section(n_elements,
                                     type,
                                     face_index,
                                     face_num,
                                     vertex_index,
                                     vertex_num,
                                     parent_element_num);

  this_nodal->sections[n_sections] = new_section;
  this_nodal->n_sections += 1;

  /* Update main structure information */

  switch(new_section->entity_dim) {
  case 3:
    this_nodal->n_cells += n_elements;
    break;
  case 2:
    this_nodal->n_faces += n_elements;
    break;
  case 1:
    this_nodal->n_edges += n_elements;
    break;
  default:
    assert(0);
  }

}

/*----------------------------------------------------------------------------
 * Append a new section to an existing fvm_nodal mesh, sharing the given
 * given connectivity and optional parent number arrays with the caller.
 *
 * The caller should not destroy or modify the arrays passed to this
 * function until the nodal mesh is destroyed.
 *
 * parameters:
 *   this_nodal         <-> nodal mesh structure
 *   n_elements         <-- number of elements to add
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   parent_element_num <-- element -> parent element number (1 to n) if non
 *                          trivial (i.e. if element definitions correspond
 *                          to a subset of the parent mesh), NULL otherwise
 *----------------------------------------------------------------------------*/

void
fvm_nodal_append_shared(fvm_nodal_t    *this_nodal,
                        cs_lnum_t       n_elements,
                        fvm_element_t   type,
                        cs_lnum_t       face_index[],
                        cs_lnum_t       face_num[],
                        cs_lnum_t       vertex_index[],
                        cs_lnum_t       vertex_num[],
                        cs_lnum_t       parent_element_num[])
{
  fvm_nodal_section_t  *new_section = NULL;
  int  n_sections = 0;

  assert(this_nodal != NULL);

  n_sections = this_nodal->n_sections;

  /* Create new section */

  BFT_REALLOC(this_nodal->sections, n_sections + 1, fvm_nodal_section_t *);

  new_section = _map_to_section(n_elements,
                                type,
                                face_index,
                                face_num,
                                vertex_index,
                                vertex_num,
                                parent_element_num);

  this_nodal->sections[n_sections] = new_section;
  this_nodal->n_sections += 1;

  /* Update main structure information */

  switch(new_section->entity_dim) {
  case 3:
    this_nodal->n_cells += n_elements;
    break;
  case 2:
    this_nodal->n_faces += n_elements;
    break;
  case 1:
    this_nodal->n_edges += n_elements;
    break;
  default:
    assert(0);
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
