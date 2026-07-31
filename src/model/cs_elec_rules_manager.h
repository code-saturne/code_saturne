#ifndef CS_ELEC_RULES_MANAGER_H
#define CS_ELEC_RULES_MANAGER_H

/*============================================================================
 * Class to parse and manage ElectricalRules.xml
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
#include <vector>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_tree.h"

/*============================================================================
 * Structure definitions
 *============================================================================*/

/* Joule model entry */
struct cs_elec_joule_model_t {
  std::string name;         /* "AC/DC", "three-phase"... */
  int         ieljou;       /* -1, 1, 2, 3, 4 */
  std::string description;
};

/* Electric arc model entry */
struct cs_elec_arc_model_t {
  std::string name;         /* "off", "on" */
  int         ielarc;       /* -1, 2 */
  std::string description;
};

/* Physical property entry */
struct cs_elec_property_t {
  std::string name;
  std::string label;
  std::string unit;
  std::string model;  /* "arc", "joule_transformer", or "" for all */
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_elec_rules_manager {

private:
  cs_tree_node_t *rules_tree_;

  std::map<std::string, cs_elec_joule_model_t>  joule_models_;
  std::map<std::string, cs_elec_arc_model_t>    arc_models_;
  std::map<std::string, std::string>            defaults_;
  std::vector<cs_elec_property_t>               properties_;
  std::vector<int>                              allowed_ieljou_;
  std::vector<int>                              allowed_ielarc_;

  void parse_rules_();

public:
  cs_elec_rules_manager(const char  *rules_xml_path);
  ~cs_elec_rules_manager();

  /* Get ieljou value from joule model name */
  int
  get_ieljou(const std::string  &name) const;

  /* Get ielarc value from arc model name */
  int
  get_ielarc(const std::string  &name) const;

  /* Check if ieljou value is valid */
  bool
  is_valid_ieljou(int  ieljou) const;

  /* Check if ielarc value is valid */
  bool
  is_valid_ielarc(int  ielarc) const;

  /* Get default value */
  std::string
  get_default(const std::string  &key) const;
  double
  get_default_double(const std::string  &key) const;

  /* Get properties for a given model */
  std::vector<cs_elec_property_t>
  get_properties(const std::string  &model) const;
};

/*============================================================================
 * Singleton access
 *============================================================================*/

cs_elec_rules_manager *
cs_get_elec_rules_manager(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_ELEC_RULES_MANAGER_H */
