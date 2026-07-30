#ifndef CS_BOUNDARY_CONDITIONS_RULES_MANAGER_H
#define CS_BOUNDARY_CONDITIONS_RULES_MANAGER_H

/*============================================================================
 * Class to parse and manage BoundaryConditionsRules.xml
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

#include <map>
#include <string>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_tree.h"

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* BC code mapping: boundary_type -> icodcl */
struct cs_bc_code_entry_t {
  std::string boundary_type;  /* e.g. "CS_BOUNDARY_WALL" */
  int         icodcl;         /* legacy solver code: 1, 5, 6, 13... */
  std::string label;          /* human-readable label */
};

/* Legacy type mapping: nature -> legacy code */
struct cs_bc_legacy_entry_t {
  std::string nature;       /* "wall", "symmetry", "default" */
  std::string legacy_code;  /* "CS_SMOOTHWALL", "CS_SYMMETRY" */
  int         value;        /* numeric value */
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_boundary_conditions_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  std::map<std::string, cs_bc_code_entry_t>  bc_code_map_;
  std::map<std::string, cs_bc_legacy_entry_t> legacy_map_;
  int default_icodcl_;

  void parse_rules_();

public:
  cs_boundary_conditions_rules_manager(const char *rules_xml_path);
  ~cs_boundary_conditions_rules_manager();

  /* Get icodcl for a given boundary type flag name */
  int
  get_icodcl(const std::string  &boundary_type) const;

  /* Get legacy type code for a given nature */
  int
  get_legacy_type(const std::string  &nature) const;
};

/*============================================================================
 * Singleton access
 *============================================================================*/

cs_boundary_conditions_rules_manager *
cs_get_boundary_conditions_rules_manager(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_BOUNDARY_CONDITIONS_RULES_MANAGER_H */
