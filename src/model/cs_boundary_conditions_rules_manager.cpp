/*============================================================================
 * Implementation of the parser for BoundaryConditionsRules.xml
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

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "base/cs_base.h"
#include "gui/cs_tree_xml.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions_rules_manager.h"

/*============================================================================
 * Private methods
 *============================================================================*/

void
cs_boundary_conditions_rules_manager::parse_rules_()
{
  /* Parse BCCodeMappings */
  cs_tree_node_t *bcm = cs_tree_node_get_child(rules_tree_,
                                                "BCCodeMappings");
  if (bcm == nullptr) return;

  /* Parse BCCodeMapping */
  cs_tree_node_t *mapping = cs_tree_node_get_child(bcm, "BCCodeMapping");
  if (mapping != nullptr) {
    for (cs_tree_node_t *entry = cs_tree_node_get_child(mapping, "Entry");
         entry != nullptr;
         entry = cs_tree_node_get_next_of_name(entry)) {

      const char *bt     = cs_tree_node_get_tag(entry, "BoundaryType");
      const char *icodcl = cs_tree_node_get_tag(entry, "icodcl");
      const char *label  = cs_tree_node_get_tag(entry, "Label");

      if (bt == nullptr || icodcl == nullptr) continue;

      cs_bc_code_entry_t e;
      e.boundary_type = bt;
      e.icodcl        = atoi(icodcl);
      e.label         = label ? label : "";

      if (e.boundary_type == "default")
        default_icodcl_ = e.icodcl;
      else
        bc_code_map_[e.boundary_type] = e;
    }
  }

  /* Parse LegacyTypeMapping */
  cs_tree_node_t *ltm = cs_tree_node_get_child(bcm, "LegacyTypeMapping");
  if (ltm != nullptr) {
    for (cs_tree_node_t *entry = cs_tree_node_get_child(ltm, "Entry");
         entry != nullptr;
         entry = cs_tree_node_get_next_of_name(entry)) {

      const char *nature      = cs_tree_node_get_tag(entry, "Nature");
      const char *legacy_code = cs_tree_node_get_tag(entry, "LegacyCode");
      const char *value       = cs_tree_node_get_tag(entry, "Value");

      if (nature == nullptr || value == nullptr) continue;

      cs_bc_legacy_entry_t e;
      e.nature      = nature;
      e.legacy_code = legacy_code ? legacy_code : "";
      e.value       = atoi(value);

      legacy_map_[e.nature] = e;
    }
  }
}

/*============================================================================
 * Constructor / Destructor
 *============================================================================*/

cs_boundary_conditions_rules_manager::cs_boundary_conditions_rules_manager
(
  const char  *rules_xml_path
)
  : default_icodcl_(1)
{
  rules_tree_ = cs_tree_node_create("");
  cs_tree_xml_read(rules_tree_, rules_xml_path);
  if (rules_tree_ == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Cannot load BoundaryConditionsRules.xml: %s",
              rules_xml_path);
  parse_rules_();
}

cs_boundary_conditions_rules_manager::~cs_boundary_conditions_rules_manager()
{
  if (rules_tree_ != nullptr)
    cs_tree_node_free(&rules_tree_);
}

/*============================================================================
 * Public methods
 *============================================================================*/

int
cs_boundary_conditions_rules_manager::get_icodcl
(
  const std::string  &boundary_type
) const
{
  auto it = bc_code_map_.find(boundary_type);
  if (it != bc_code_map_.end())
    return it->second.icodcl;
  return default_icodcl_;
}

int
cs_boundary_conditions_rules_manager::get_legacy_type
(
  const std::string  &nature
) const
{
  auto it = legacy_map_.find(nature);
  if (it != legacy_map_.end())
    return it->second.value;
  /* default */
  auto def = legacy_map_.find("default");
  if (def != legacy_map_.end())
    return def->second.value;
  return 0;
}

/*============================================================================
 * Singleton
 *============================================================================*/

cs_boundary_conditions_rules_manager *
cs_get_boundary_conditions_rules_manager(bool  no_instanciate)
{
  static cs_boundary_conditions_rules_manager *instance = nullptr;
  if (instance == nullptr && no_instanciate == false) {
    const char *datadir = cs_base_get_pkgdatadir();
    char path[1024];
    snprintf(path, sizeof(path) - 1,
             "%s/model/BoundaryConditionsRules.xml", datadir);
    path[sizeof(path) - 1] = '\0';
    instance = new cs_boundary_conditions_rules_manager(path);
  }
  return instance;
}

/*----------------------------------------------------------------------------*/
