/*============================================================================
 * Parse TurbulentFluxRules.xml and expose rules for creation of
 * turbulente flux fields.
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
#include <stdexcept>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "bft/bft_error.h"
#include "gui/cs_tree_xml.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulent_flux_rules_manager.h"

/*============================================================================
 * Static global variable (singleton instance)
 *============================================================================*/

static cs_turbulent_flux_rules_manager *_instance = nullptr;

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
  if (_instance != nullptr) {
    delete _instance;
    _instance = nullptr;
  }
}

/*============================================================================
 * Static helper methods
 *============================================================================*/

inline int
cs_turbulent_flux_rules_manager::_atoi_safe(const char *s, int def)
{
  if (s == nullptr || s[0] == '\0') return def;
  return atoi(s);
}

/*============================================================================
 * Private methods
 *============================================================================*/

cs_turb_flux_field_t
cs_turbulent_flux_rules_manager::parse_field_(cs_tree_node_t *fn)
{
  cs_turb_flux_field_t f;

  const char *name_t  = cs_tree_node_get_tag(fn, "NameTemplate");
  const char *label_t = cs_tree_node_get_tag(fn, "LabelTemplate");
  const char *dim     = cs_tree_node_get_tag(fn, "Dimension");
  const char *type    = cs_tree_node_get_tag(fn, "Type");
  const char *coupled = cs_tree_node_get_tag(fn, "Coupled");
  const char *clip    = cs_tree_node_get_tag(fn, "Clipping");
  const char *post    = cs_tree_node_get_tag(fn, "PostProcess");
  const char *log     = cs_tree_node_get_tag(fn, "Log");

  f.name_template  = name_t  ? name_t  : "";
  f.label_template = label_t ? label_t : "";
  f.dimension      = _atoi_safe(dim, 1);
  f.type           = type    ? type    : "property";
  f.coupled        = coupled && atoi(coupled) == 1;
  f.clipping       = clip    && atoi(clip)    == 1;
  f.post_process   = post    && atoi(post)    == 1;
  f.log            = log     && atoi(log)     == 1;

  return f;
}

void
cs_turbulent_flux_rules_manager::parse_rules_()
{
  cs_tree_node_t *fcr = cs_tree_node_get_child(rules_tree_,
                                               "FieldCreationRules");
  if (fcr == nullptr) return;

  for (cs_tree_node_t *rule = cs_tree_node_get_child(fcr, "Rule");
       rule != nullptr;
       rule = cs_tree_node_get_next_of_name(rule)) {

    const char *model = cs_tree_node_get_tag(rule, "Model");
    const char *num   = cs_tree_node_get_tag(rule, "NumericValue");

    if (model == nullptr) continue;

    cs_turb_flux_rule_t r;
    r.model_name    = model;
    r.numeric_value = _atoi_safe(num, 0);

    const char *desc = cs_tree_node_get_child_value_str(rule, "Description");
    if (desc) r.description = desc;

    for (cs_tree_node_t *fn = cs_tree_node_get_child(rule, "Field");
         fn != nullptr;
         fn = cs_tree_node_get_next_of_name(fn)) {
      r.fields.push_back(parse_field_(fn));
    }

    rules_by_model_[r.model_name]    = r;
    rules_by_value_[r.numeric_value] = r;
  }
}

/*============================================================================
 * Constructor / Destructor
 *============================================================================*/

cs_turbulent_flux_rules_manager::cs_turbulent_flux_rules_manager
  (const char *rules_xml_path)
{
  rules_tree_ = cs_tree_node_create("");
  cs_tree_xml_read(rules_tree_, rules_xml_path);
  if (rules_tree_ == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "Cannot load TurbulentFluxRules.xml: %s", rules_xml_path);

  parse_rules_();
}

cs_turbulent_flux_rules_manager::~cs_turbulent_flux_rules_manager()
{
  if (rules_tree_ != nullptr)
    cs_tree_node_free(&rules_tree_);
}

/*============================================================================
 * Public methods
 *============================================================================*/

const cs_turb_flux_rule_t *
cs_turbulent_flux_rules_manager::get_rule_by_model
  (const std::string &model) const
{
  auto it = rules_by_model_.find(model);
  if (it != rules_by_model_.end())
    return &it->second;
  return nullptr;
}

const cs_turb_flux_rule_t *
cs_turbulent_flux_rules_manager::get_rule_by_value(int numeric_value) const
{
  auto it = rules_by_value_.find(numeric_value);
  if (it != rules_by_value_.end())
    return &it->second;
  return nullptr;
}

std::string
cs_turbulent_flux_rules_manager::resolve_name
  (const std::string &name_template, const std::string &scalar_name)
{
  std::string result = name_template;
  size_t pos = result.find("{scalar}");
  while (pos != std::string::npos) {
    result.replace(pos, 8, scalar_name);
    pos = result.find("{scalar}", pos + scalar_name.size());
  }
  return result;
}

int
cs_turbulent_flux_rules_manager::get_numeric_value
  (const std::string &model_name) const
{
  auto it = rules_by_model_.find(model_name);
  if (it != rules_by_model_.end())
    return it->second.numeric_value;
  return 0; /* SGDH par defaut */
}

bool
cs_turbulent_flux_rules_manager::creates_fields
  (const std::string &model_name) const
{
  auto it = rules_by_model_.find(model_name);
  if (it == rules_by_model_.end()) return false;
  return !it->second.fields.empty();
}

bool
cs_turbulent_flux_rules_manager::creates_alpha_field
  (const std::string &model_name) const
{
  auto it = rules_by_model_.find(model_name);
  if (it == rules_by_model_.end()) return false;
  for (const auto &f : it->second.fields) {
    if (f.name_template.find("_alpha") != std::string::npos)
      return true;
  }
  return false;
}

/*============================================================================
 * Singleton
 *============================================================================*/

cs_turbulent_flux_rules_manager *
cs_get_turbulent_flux_rules_manager(void)
{
  if (_instance == nullptr) {
    const char *datadir = cs_base_get_pkgdatadir();
    char path[1024];
    snprintf(path, sizeof(path) - 1,
             "%s/model/TurbulentFluxRules.xml", datadir);
    path[sizeof(path) - 1] = '\0';
    _instance = new cs_turbulent_flux_rules_manager(path);

    cs_base_at_finalize(_rule_manager_finalize);
  }

  return _instance;
}

/*-----------------------------------------------------------------------------*/
