/*============================================================================
 * Implementation of the parser for LagrangianRules.xml
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

#include "cs_lagr_rules_manager.h"

/*============================================================================
 * Private methods
 *============================================================================*/

void
cs_lagr_rules_manager::parse_rules_()
{
  cs_tree_node_t *defs = cs_tree_node_get_child(rules_tree_, "Definitions");
  if (defs == nullptr) return;

  /* Parse coupling modes */
  for (cs_tree_node_t *td = cs_tree_node_get_child(defs, "TypeDef");
       td != nullptr;
       td = cs_tree_node_get_next_of_name(td)) {

    const char *tname = cs_tree_node_get_tag(td, "name");
    if (tname == nullptr) continue;

    if (std::string(tname) == "LagrangianCouplingMode") {
      for (cs_tree_node_t *v = cs_tree_node_get_child(td, "Value");
           v != nullptr;
           v = cs_tree_node_get_next_of_name(v)) {
        const char *name  = cs_tree_node_get_tag(v, "name");
        const char *en    = cs_tree_node_get_tag(v, "Enum");
        const char *iv    = cs_tree_node_get_tag(v, "Int");
        const char *desc  = cs_tree_node_get_tag(v, "Description");
        if (name == nullptr) continue;
        cs_lagr_coupling_mode_t m;
        m.name        = name;
        m.enum_name   = en   ? en   : "";
        m.int_value   = iv   ? atoi(iv) : 0;
        m.description = desc ? desc : "";
        coupling_modes_[m.name] = m;
      }
    }
    else if (std::string(tname) == "ParticlesModel") {
      for (cs_tree_node_t *v = cs_tree_node_get_child(td, "Value");
           v != nullptr;
           v = cs_tree_node_get_next_of_name(v)) {
        const char *name  = cs_tree_node_get_tag(v, "name");
        const char *en    = cs_tree_node_get_tag(v, "Enum");
        const char *iv    = cs_tree_node_get_tag(v, "Int");
        const char *desc  = cs_tree_node_get_tag(v, "Description");
        if (name == nullptr) continue;
        cs_lagr_phys_model_t m;
        m.name        = name;
        m.enum_name   = en   ? en   : "";
        m.int_value   = iv   ? atoi(iv) : 0;
        m.description = desc ? desc : "";
        phys_models_[m.name] = m;
      }
    }
  }

  /* Parse defaults */
  cs_tree_node_t *dflts = cs_tree_node_get_child(rules_tree_, "Defaults");
  if (dflts != nullptr) {
    for (cs_tree_node_t *d = cs_tree_node_get_child(dflts, "Default");
         d != nullptr;
         d = cs_tree_node_get_next_of_name(d)) {
      const char *key = cs_tree_node_get_tag(d, "Key");
      const char *val = cs_tree_node_get_tag(d, "Value");
      if (key != nullptr && val != nullptr)
        defaults_[key] = val;
    }
  }
}

/*============================================================================
 * Constructor / Destructor
 *============================================================================*/

cs_lagr_rules_manager::cs_lagr_rules_manager(const char *rules_xml_path)
{
  rules_tree_ = cs_tree_node_create("");
  cs_tree_xml_read(rules_tree_, rules_xml_path);
  if (rules_tree_ == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Cannot load LagrangianRules.xml: %s", rules_xml_path);
  parse_rules_();
}

cs_lagr_rules_manager::~cs_lagr_rules_manager()
{
  if (rules_tree_ != nullptr)
    cs_tree_node_free(&rules_tree_);
}

/*============================================================================
 * Public methods
 *============================================================================*/

int
cs_lagr_rules_manager::get_coupling_mode(const std::string &name) const
{
  auto it = coupling_modes_.find(name);
  if (it != coupling_modes_.end())
    return it->second.int_value;
  return 0;
}

int
cs_lagr_rules_manager::get_phys_model(const std::string &name) const
{
  auto it = phys_models_.find(name);
  if (it != phys_models_.end())
    return it->second.int_value;
  return 0;
}

std::string
cs_lagr_rules_manager::get_default(const std::string &key) const
{
  auto it = defaults_.find(key);
  if (it != defaults_.end())
    return it->second;
  return "";
}

double
cs_lagr_rules_manager::get_default_double(const std::string &key) const
{
  std::string val = get_default(key);
  if (!val.empty())
    return atof(val.c_str());
  return 0.0;
}

int
cs_lagr_rules_manager::get_default_int(const std::string &key) const
{
  std::string val = get_default(key);
  if (!val.empty())
    return atoi(val.c_str());
  return 0;
}

bool
cs_lagr_rules_manager::requires_restart(const std::string &coupling_mode) const
{
  return coupling_mode == "frozen";
}

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_lagr_rules_manager *_g_lagr_rules_manager = nullptr;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy helper and associated tree.
 */
/*----------------------------------------------------------------------------*/

static void
_rule_manager_finalize(void)
{
  if (_g_lagr_rules_manager != nullptr) {
    delete _g_lagr_rules_manager;
    _g_lagr_rules_manager = nullptr;
  }
}

/*============================================================================
 * Singleton
 *============================================================================*/

cs_lagr_rules_manager *
cs_get_lagr_rules_manager(void)
{
  if (_g_lagr_rules_manager == nullptr) {
    const char *datadir = cs_base_get_pkgdatadir();
    char path[1024];
    snprintf(path, sizeof(path) - 1,
             "%s/model/LagrangianRules.xml", datadir);
    path[sizeof(path) - 1] = '\0';
    _g_lagr_rules_manager = new cs_lagr_rules_manager(path);
    cs_base_at_finalize(_rule_manager_finalize);
  }
  return _g_lagr_rules_manager;
}

/*----------------------------------------------------------------------------*/
