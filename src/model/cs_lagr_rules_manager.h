#ifndef CS_LAGR_RULES_MANAGER_H
#define CS_LAGR_RULES_MANAGER_H

/*============================================================================
 * Class to parse and manage LagrangianRules.xml
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

struct cs_lagr_coupling_mode_t {
  std::string name;         /* "off", "one_way", "two_way", "frozen" */
  std::string enum_name;    /* "CS_LAGR_OFF", "CS_LAGR_ONEWAY_COUPLING"... */
  int         int_value;    /* 0, 1, 2, 3 */
  std::string description;
};

struct cs_lagr_phys_model_t {
  std::string name;         /* "off", "thermal", "coal", "ctwr" */
  std::string enum_name;    /* "CS_LAGR_PHYS_OFF"... */
  int         int_value;    /* 0, 1, 2, 3 */
  std::string description;
};

/*============================================================================
 * Class definition
 *============================================================================*/

class cs_lagr_rules_manager {
private:
  cs_tree_node_t *rules_tree_;

  std::map<std::string, cs_lagr_coupling_mode_t> coupling_modes_;
  std::map<std::string, cs_lagr_phys_model_t>    phys_models_;
  std::map<std::string, std::string>              defaults_;

  void parse_rules_();

public:
  cs_lagr_rules_manager(const char *rules_xml_path);
  ~cs_lagr_rules_manager();

  /* Get coupling mode int value from name */
  int get_coupling_mode(const std::string &name) const;

  /* Get physical model int value from name */
  int get_phys_model(const std::string &name) const;

  /* Get default value */
  std::string get_default(const std::string &key) const;
  double      get_default_double(const std::string &key) const;
  int         get_default_int(const std::string &key) const;

  /* Check if model requires restart */
  bool requires_restart(const std::string &coupling_mode) const;
};

/*============================================================================
 * Singleton access
 *============================================================================*/

cs_lagr_rules_manager *cs_get_lagr_rules_manager(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_LAGR_RULES_MANAGER_H */
