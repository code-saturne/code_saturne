/*============================================================================
 * cs_elec_rules_manager.cpp
 *
 * Implementation of the parser for ElectricalRules.xml
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

#include "cs_elec_rules_manager.h"

/*============================================================================
 * Private methods
 *============================================================================*/

void
cs_elec_rules_manager::parse_rules_()
{
  cs_tree_node_t *defs = cs_tree_node_get_child(rules_tree_, "Definitions");
  if (defs == nullptr) return;

  /* Parse JouleModel */
  for (cs_tree_node_t *td = cs_tree_node_get_child(defs, "TypeDef");
       td != nullptr;
       td = cs_tree_node_get_next_of_name(td)) {

    const char *tname = cs_tree_node_get_tag(td, "name");
    if (tname == nullptr) continue;

    if (std::string(tname) == "JouleModel") {
      for (cs_tree_node_t *v = cs_tree_node_get_child(td, "Value");
           v != nullptr;
           v = cs_tree_node_get_next_of_name(v)) {
        const char *name = cs_tree_node_get_tag(v, "name");
        const char *iv   = cs_tree_node_get_tag(v, "Int");
        const char *desc = cs_tree_node_get_tag(v, "Description");
        if (name == nullptr) continue;
        cs_elec_joule_model_t m;
        m.name        = name;
        m.ieljou      = iv ? atoi(iv) : -1;
        m.description = desc ? desc : "";
        joule_models_[m.name] = m;
        allowed_ieljou_.push_back(m.ieljou);
      }
    }
    else if (std::string(tname) == "ElectricArcModel") {
      for (cs_tree_node_t *v = cs_tree_node_get_child(td, "Value");
           v != nullptr;
           v = cs_tree_node_get_next_of_name(v)) {
        const char *name = cs_tree_node_get_tag(v, "name");
        const char *iv   = cs_tree_node_get_tag(v, "Int");
        const char *desc = cs_tree_node_get_tag(v, "Description");
        if (name == nullptr) continue;
        cs_elec_arc_model_t m;
        m.name        = name;
        m.ielarc      = iv ? atoi(iv) : -1;
        m.description = desc ? desc : "";
        arc_models_[m.name] = m;
        allowed_ielarc_.push_back(m.ielarc);
      }
    }
  }

  /* Parse Defaults */
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

  /* Parse PhysicalProperties */
  cs_tree_node_t *props = cs_tree_node_get_child(rules_tree_,
                                                  "PhysicalProperties");
  if (props != nullptr) {
    for (cs_tree_node_t *p = cs_tree_node_get_child(props, "Property");
         p != nullptr;
         p = cs_tree_node_get_next_of_name(p)) {
      const char *name  = cs_tree_node_get_tag(p, "name");
      const char *label = cs_tree_node_get_tag(p, "Label");
      const char *unit  = cs_tree_node_get_tag(p, "Unit");
      const char *model = cs_tree_node_get_tag(p, "Model");
      if (name == nullptr) continue;
      cs_elec_property_t prop;
      prop.name  = name;
      prop.label = label ? label : "";
      prop.unit  = unit  ? unit  : "";
      prop.model = model ? model : "";
      properties_.push_back(prop);
    }
  }
}

/*============================================================================
 * Constructor / Destructor
 *============================================================================*/

cs_elec_rules_manager::cs_elec_rules_manager
(
  const char  *rules_xml_path
)
{
  rules_tree_ = cs_tree_node_create("");
  cs_tree_xml_read(rules_tree_, rules_xml_path);
  if (rules_tree_ == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Cannot load ElectricalRules.xml: %s", rules_xml_path);
  parse_rules_();
}

cs_elec_rules_manager::~cs_elec_rules_manager()
{
  if (rules_tree_ != nullptr)
    cs_tree_node_free(&rules_tree_);
}

/*============================================================================
 * Public methods
 *============================================================================*/

int
cs_elec_rules_manager::get_ieljou(const std::string &name) const
{
  auto it = joule_models_.find(name);
  if (it != joule_models_.end())
    return it->second.ieljou;
  return -1;
}

int
cs_elec_rules_manager::get_ielarc(const std::string &name) const
{
  auto it = arc_models_.find(name);
  if (it != arc_models_.end())
    return it->second.ielarc;
  return -1;
}

bool
cs_elec_rules_manager::is_valid_ieljou(int ieljou) const
{
  for (int v : allowed_ieljou_)
    if (v == ieljou) return true;
  return false;
}

bool
cs_elec_rules_manager::is_valid_ielarc(int ielarc) const
{
  for (int v : allowed_ielarc_)
    if (v == ielarc) return true;
  return false;
}

std::string
cs_elec_rules_manager::get_default(const std::string &key) const
{
  auto it = defaults_.find(key);
  if (it != defaults_.end())
    return it->second;
  return "";
}

double
cs_elec_rules_manager::get_default_double(const std::string &key) const
{
  std::string val = get_default(key);
  if (!val.empty())
    return atof(val.c_str());
  return 0.0;
}

std::vector<cs_elec_property_t>
cs_elec_rules_manager::get_properties(const std::string &model) const
{
  std::vector<cs_elec_property_t> result;
  for (const auto &p : properties_) {
    if (p.model.empty() || p.model == model)
      result.push_back(p);
  }
  return result;
}

/*============================================================================
 * Singleton
 *============================================================================*/

cs_elec_rules_manager *
cs_get_elec_rules_manager(void)
{
  static cs_elec_rules_manager *instance = nullptr;
  if (instance == nullptr) {
    const char *datadir = cs_base_get_pkgdatadir();
    char path[1024];
    snprintf(path, sizeof(path) - 1,
             "%s/model/ElectricalRules.xml", datadir);
    path[sizeof(path) - 1] = '\0';
    instance = new cs_elec_rules_manager(path);
  }
  return instance;
}

/*----------------------------------------------------------------------------*/
